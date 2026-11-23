#include "scriptAnalysis.hpp"

#include <sstream>

std::optional<SourceRange> scriptCompletionRange(const ParsedScript& script,
                                                const size_t byteOffset)
{
  std::optional<SourceRange> result;
  const auto text=std::string_view(script.source.text);
  bool applicable=byteOffset<=text.size();
  for(const auto& error : script.lexicalErrors)
    if(error.begin.offset<=byteOffset and byteOffset<=error.end.offset)
      applicable=false;
  for(const auto& trivia : script.trivia)
    {
      const auto spelling=text.substr(trivia.begin.offset,trivia.end.offset-trivia.begin.offset);
      if(trivia.begin.offset<=byteOffset and
         (byteOffset<trivia.end.offset or (byteOffset==trivia.end.offset and spelling.starts_with("//"))) and
         spelling.find_first_not_of(" \t\r\n\f\v")!=std::string_view::npos)
        applicable=false;
    }
  if(applicable)
    {
      SourceCursor cursor;
      for(const auto& token : script.tokens)
        if(not result and token.range.end.offset>=byteOffset)
          {
            const bool within=token.range.begin.offset<=byteOffset;
            if(within and token.identifier)
              result=token.range;
            else if(within and not token.end and
                    (byteOffset<token.range.end.offset or
                     (token.range.begin.offset<text.size() and text[token.range.begin.offset]>='0' and
                      text[token.range.begin.offset]<='9')))
              applicable=false;
            else
              {
                for(size_t i=cursor.offset;i<byteOffset;i++)
                  cursor.advance(text[i]);
                result=SourceRange{cursor,cursor};
              }
          }
        else if(token.range.end.offset<byteOffset)
          cursor=token.range.end;
    }
  if(not applicable)
    result.reset();
  
  return result;
}

ParsedScript parseScriptForEditor(const ScriptSource& source,
                                  const std::optional<size_t> completionOffset)
{
  static const auto sharedGrammar=std::make_shared<const pp::Grammar>(getParser());
  ParsedScript result{.source=source,.grammar=sharedGrammar};
  const auto& grammar=*result.grammar;
  const auto locate=[&result](const SourceCursor& cursor)
  {
    return result.source.locate(cursor);
  };
  
  // Recovery starts in the real lexer: invalid text becomes a diagnostic and
  // tokenization continues, preserving the positions of the remaining tokens
  auto parsed=tokenizeInput(grammar,result.source.text,locate,true);
  
  // NICE supplies possible inserted tokens, while ParsePact decides whether
  // they are legal using the same LR tables as strict parsing
  // Resolve spellings through the lexer so no terminal numbers are hardcoded
  const auto identifier=grammar.regexMatcher.match("niceEditorIdentifier");
  const size_t identifierSymbol=grammar.iSymbolOfRegex[identifier->iToken];
  std::vector<ParseToken> repairs;
  
  // Prefer ending a statement or closing an unfinished construct before
  // supplying missing syntax inside it, e.g. x=1 y=2; can insert ; before y
  // Candidate order breaks ties between repairs with the same insertion count
  for(const std::string_view spelling : {";", ")", "]", "}", "=", "{"})
    if(const auto token=grammar.regexMatcher.match(spelling);token and token->matchedString==spelling)
      repairs.push_back({grammar.iSymbolOfRegex[token->iToken],ShiftNode{.txt=spelling,.range=SourceRange{}}});
  
  // A missing identifier can stand in for an unfinished name or expression
  // Its private spelling distinguishes recovery scaffolding from user names
  repairs.push_back({identifierSymbol,ShiftNode{.txt=recoveryMarker,.range=SourceRange{}}});
  
  for(const auto& symbol : grammar.symbols)
    if(symbol.type==BaseGrammarSymbol::Type::TERMINAL_SYMBOL)
      {
        std::string_view name=symbol.name;
        while(not name.empty())
          {
            const size_t separator=name.find('|');
            const auto alternative=name.substr(0,separator);
            if(not alternative.empty() and std::all_of(alternative.begin(),alternative.end(),[](const char c)
               {
                 return (c>='a' and c<='z') or (c>='A' and c<='Z') or c=='_';
               }))
              result.keywords.emplace_back(alternative);
            name=separator==std::string_view::npos?std::string_view{}:name.substr(separator+1);
          }
      }
  
  for(const auto& token : parsed.tokens)
    result.tokens.push_back({*token.shift.range,token.symbol==identifierSymbol,token.symbol==grammar.iEndSymbol});
  for(const auto& trivia : parsed.trivia)
    result.trivia.push_back(*trivia.range);
  for(const auto& error : parsed.diagnostics)
    result.lexicalErrors.push_back(error.range);
  
  if(completionOffset)
    if(const auto range=scriptCompletionRange(result,*completionOffset))
      {
        result.completion=range;
        size_t index=0;
        while(index<parsed.tokens.size() and
              (parsed.tokens[index].shift.range->end.offset<range->begin.offset or
               (range->begin!=range->end and parsed.tokens[index].shift.range->end.offset==range->begin.offset)))
          index++;
        if(range->begin!=range->end)
          parsed.tokens[index].shift.txt=completionMarker;
        else
          {
            if(index<parsed.tokens.size() and parsed.tokens[index].symbol!=grammar.iEndSymbol and
               parsed.tokens[index].shift.range->end.offset==range->begin.offset)
              index++;
            const auto cursor=range->begin;
            parsed.tokens.insert(parsed.tokens.begin()+index,
              {identifierSymbol,ShiftNode{.txt=completionMarker,.cursor=cursor,.location=locate(cursor),.range=*range}});
          }
      }
  
  // At an LR error ParsePact records the offending token, then forecasts parser
  // states without building trees or executing semantic actions
  // A breadth-first search tries up to four insertions and 128 plans, choosing
  // the first sequence that also lets the original token be consumed
  // Successful insertions are replayed into the real parser at zero-width
  // ranges before that token; if none works, the offending token is discarded
  // At EOF an individually legal insertion may instead advance the parser
  // toward completion over several rounds; all recovery stops after 256 rounds
  // These limits belong to ParsePact's parseTokens implementation
  //
  // Repair is a best-effort structural guess, not a correction of the script's
  // meaning: diagnostics remain visible and later errors can be consequences
  // of an earlier guess, but usable declarations can survive on either side
  parseTokens(grammar,parsed,repairs,true);
  
  // Convert the completed tree or surviving forest with the normal NICE AST
  // actions in editor mode, recording node ranges for the completion visitor
  // Constructing this AST does not evaluate script or invoke registered hosts
  const auto executor=getParseTreeExecuctor(grammar.actionStrings,true);
  for(const auto& tree : parsed.roots)
    result.roots.push_back(executor.execParseTree(flattenParseTree(tree),
      [&result](const std::shared_ptr<ASTNode>& node,
                const SourceRange& range)
      {
        result.ranges[node.get()]=range;
      }));
  result.diagnostics=std::move(parsed.diagnostics);
  std::stable_sort(result.diagnostics.begin(),result.diagnostics.end(),[](const auto& left,
                                                                         const auto& right)
  {
    return left.range.begin.offset<right.range.begin.offset;
  });
  result.diagnostics.erase(std::unique(result.diagnostics.begin(),result.diagnostics.end(),[](const auto& left,
                                                                                            const auto& right)
  {
    return left.range.begin==right.range.begin and left.range.end==right.range.end and left.message==right.message;
  }),result.diagnostics.end());
  
  return result;
}

void writeScriptLispString(std::ostream& output,
                           const std::string_view value)
{
  output<<'"';
  for(const unsigned char c : value)
    {
      if(c=='"' or c=='\\')
        output<<'\\'<<c;
      else if(c<32 or c==127)
        output<<'\\'<<char('0'+((c>>6)&7))<<char('0'+((c>>3)&7))<<char('0'+(c&7));
      else
        output<<c;
    }
  output<<'"';
}

std::string formatScriptDiagnosticsForEmacs(const ParsedScript& script)
{
  std::ostringstream output;
  output<<'(';
  for(const auto& diagnostic : script.diagnostics)
    {
      const auto begin=script.source.locate(diagnostic.range.begin);
      const auto end=script.source.locate(diagnostic.range.end);
      output<<'(';
      writeScriptLispString(output,begin.source);
      output<<' '<<begin.cursor.offset<<' '<<std::max(begin.cursor.offset,end.cursor.offset)
            <<' '<<begin.cursor.line<<' '<<begin.cursor.column<<' ';
      writeScriptLispString(output,diagnostic.message);
      output<<')';
    }
  output<<')';
  
  return output.str();
}
