#include "nASTy.hpp"

#include <cmath>
#include <cstdlib>
#include <charconv>
#include <sstream>
#include <format>

#include "routines/ios.hpp"

std::string precedingDocumentation(const pp::internal::ShiftNode& token)
{
  StringMatcher matcher(token.leadingTrivia,token.triviaCursor);
  std::string description;
  size_t lastDocLine=0;
  bool cleanLine=token.triviaCursor.column==1;
  
  while(not matcher.ref.empty())
    {
      if(matcher.ref.front()=='\r' or matcher.ref.front()=='\n')
        {
          const size_t previousLine=matcher.cursor.line;
          matcher.advance();
          
          if(matcher.cursor.line!=previousLine)
            {
              cleanLine=true;
              
              if(lastDocLine!=0 and matcher.cursor.line>lastDocLine+1)
                {
                  description.clear();
                  lastDocLine=0;
                }
            }
        }
      else if(matcher.matchAnyCharIn(" \t\f\v"))
        {
          // Indentation preserves adjacency to the following definition
        }
      else if(matcher.ref.starts_with("//"))
        {
          const bool documented=cleanLine and
            (matcher.ref.starts_with("///") or matcher.ref.starts_with("//!"));
          const std::string_view before=matcher.ref;
          const size_t line=matcher.cursor.line;
          matcher.matchLineComment();
          
          if(documented)
            {
              if(lastDocLine!=0)
                description+='\n';
              
              std::string_view text=before.substr(3,before.size()-matcher.ref.size()-3);
              if(text.starts_with(' '))
                text.remove_prefix(1);
              
              description+=text;
              lastDocLine=line;
            }
          else
            {
              description.clear();
              lastDocLine=0;
            }
          
          cleanLine=false;
        }
      else
        {
          // Ordinary block comments and any intervening text end the documentation block
          if(not matcher.matchBlockComment())
            matcher.advance();
          
          description.clear();
          lastDocLine=0;
          cleanLine=false;
        }
    }
  
  if(not cleanLine or matcher.cursor.line!=lastDocLine+1)
    description.clear();
  
  return description;
}

namespace
{
  /// Embed the standalone help text without a runtime file dependency
  constexpr char languageHelp[]=
#include "nASTyHelp.inc"
    ;
  
  /// Describe a default without evaluating script code or accessing the calling scope
  std::string describeArgumentDefault(const ASTNode& node)
  {
    std::ostringstream text;
    
    if(const auto* literal=std::get_if<ValueNode>(&node))
      {
        if(const auto* string=std::get_if<std::string>(&literal->value.data))
          {
            text<<'"';
            for(const char c : *string)
              switch(c)
                {
                case '\n': text<<"\\n"; break;
                case '\r': text<<"\\r"; break;
                case '\t': text<<"\\t"; break;
                case '\\': text<<"\\\\"; break;
                case '"': text<<"\\\""; break;
                default: text<<c;
                }
            
            text<<'"';
          }
        else if(const auto* number=std::get_if<double>(&literal->value.data))
          text<<std::format("{}",*number);
        else
          text<<literal->value;
      }
    else if(const auto* identifier=std::get_if<IdNode>(&node))
      text<<identifier->name;
    else
      text<<"<expression>";
    
    return text.str();
  }
  
  /// Mark an argument as optional and show its unevaluated default when present
  std::string describeArgument(const std::string& name,
                               const std::shared_ptr<ASTNode>& defaultValue)
  {
    std::string argument=name;
    
    if(defaultValue)
      argument="["+argument+" = "+describeArgumentDefault(*defaultValue)+"]";
    
    return argument;
  }
  
  /// Indent documentation lines and wrap prose without splitting words or losing indentation
  void writeDocumentation(std::ostream& output,
                          const std::string& description,
                          const std::string& prefix)
  {
    std::istringstream lines(description);
    std::string line;
    
    while(std::getline(lines,line))
      {
        const size_t firstText=line.find_first_not_of(" \t");
        const std::string indentation=line.substr(0,firstText==std::string::npos?line.size():firstText);
        const size_t width=88>prefix.size()?88-prefix.size():1;
        bool wrapping=true;
        
        while(wrapping and line.size()>width)
          {
            const size_t split=line.find_last_of(" \t",width);
            wrapping=split!=std::string::npos and split>indentation.size();
            
            if(wrapping)
              {
                output<<prefix<<line.substr(0,split)<<'\n';
                line=indentation+line.substr(split+1);
              }
          }
        
        output<<prefix<<line<<'\n';
      }
  }
  
  /// Write every supported signature and its documentation as one separated entry
  void writeFunctionEntry(std::ostream& output,
                          const std::string& name,
                          const EntityInfo& info)
  {
    for(const auto& arguments : info.argumentForms)
      {
        output<<"  "<<name<<'(';
        
        for(size_t i=0;i<arguments.size();i++)
          output<<(i==0?"":", ")<<arguments[i];
        
        output<<")\n";
      }
    
    output<<"    Kind: "<<info.kind<<'\n';
    
    if(info.location)
      output<<"    Location: "<<info.location->describe()<<'\n';
    
    if(not info.description.empty())
      {
        output<<"    Description:\n";
        writeDocumentation(output,info.description,"      ");
      }
  }
  
  /// Write a variable within the separate variable group
  void writeVariableEntry(std::ostream& output,
                          const std::string& name,
                          const EntityInfo& info)
  {
    output<<"  "<<name<<": "<<info.kind<<'\n';
    
    if(info.location)
      output<<"    Location: "<<info.location->describe()<<'\n';
    
    if(not info.description.empty())
      writeDocumentation(output,info.description,"    ");
  }
}

FunctionArgumentForms entityArgumentForms(const HostFunction& function)
{
  FunctionArgumentForms forms=function.argumentForms;
  
  if(forms.empty())
    {
      auto& arguments=forms.emplace_back();
      for(const auto& name : function.argumentNames)
        arguments.push_back(function.isConstructor?"["+name+"]":name);
      
      if(arguments.empty() and not function.isConstructor)
        arguments.push_back("<positional arguments; signature unavailable>");
    }
  
  return forms;
}

FunctionArgumentForms entityArgumentForms(const Function& function)
{
  FunctionArgumentForms forms(1);
  
  for(const auto& parameter : function.funcNode->pars.list)
    {
      const std::string prefix=parameter.isRef?"&":"";
      
      if(parameter.name==FuncParNode::VARIADIC_PAR_NAME)
        forms.front().push_back(prefix+"...");
      else
        forms.front().push_back(describeArgument(prefix+parameter.name,parameter.def));
    }
  
  return forms;
}

FunctionArgumentForms entityArgumentForms(const StructDef& definition)
{
  FunctionArgumentForms forms(1);
  
  for(const auto& member : definition.structDefNode->members.list)
    forms.front().push_back(describeArgument(member.name,member.def));
  
  return forms;
}

std::string formatEntityListing(const Environment& environment,
                                const EntitySelection selection,
                                const bool verbose)
{
  std::ostringstream output;
  
  if(not verbose)
    for(const auto& name : environment.entityNames(selection))
      output<<name<<'\n';
  else
    {
      for(const auto group : {EntitySelection::Functions,EntitySelection::Variables})
        if(selection==EntitySelection::All or selection==group)
          {
            const auto names=environment.entityNames(group);
            
            if(not names.empty())
              {
                const bool functions=group==EntitySelection::Functions;
                output<<(functions?"Functions":"Variables")<<"\n---------\n";
                
                if(functions)
                  output<<"  Arguments in [...] are optional; ... accepts zero or more arguments.\n";
                
                output<<'\n';
                for(const auto& name : names)
                  {
                    const EntityInfo info=environment.entityInfo(name);
                    
                    if(functions)
                      writeFunctionEntry(output,name,info);
                    else
                      writeVariableEntry(output,name,info);
                    
                    output<<'\n';
                  }
              }
          }
    }
  
  return output.str();
}

std::string formatHelp(const Environment& environment,
                       const std::string& keyword)
{
  std::ostringstream output;
  
  if(keyword.empty())
    output<<languageHelp;
  else
    {
      const auto names=environment.entityNames();
      if(std::binary_search(names.begin(),names.end(),keyword))
        {
          const EntityInfo info=environment.entityInfo(keyword);
          
          if(info.callable)
            writeFunctionEntry(output,keyword,info);
          else
            writeVariableEntry(output,keyword,info);
        }
      else
        output<<"No visible entity named '"<<keyword<<"'. Use listEntities() to find available names.\n";
    }
  
  return output.str();
}

/// Build the ParsePact parser for the script grammar and its semantic actions
/// The grammar text has static lifetime because parser data retains string views
pp::Grammar getParser()
{
#define BARE_WHITESPACES			\
  "( |\\n|\\r|\\t)+"
  
#define CPP_COMMENT				\
  "//[^\\n\\r]*"
  
#define C_COMMENT				\
  "/\\*([^\\*]|\\*+[^\\*/])*\\*+/"
  
#define LINE(ARG) ARG "\n"
  
  // The returned grammar retains string views into this storage
  static constexpr char cGrammar[]=
    LINE("c {")
    LINE(    "   %whitespace \"" BARE_WHITESPACES "|" CPP_COMMENT "|" C_COMMENT "\";")
    LINE(    "   %none lowerThanElse;")
    LINE(    "   %none \"else\";")
    LINE(    "   %left \"\\(\";")
    LINE(    "   %left \"\\)\";")
    LINE(    "   %none BRACKETS;")
    LINE(    "   %left \",\";")
    LINE(    "   %right \"=\" \"\\+=\" \"\\-=\" \"\\*=\" \"/=\";")
    LINE(    "   %left \"\\|\\||or\";")
    LINE(    "   %left \"&&|and\";")
    LINE(    "   %left \"==\" \"!=\";")
    LINE(    "   %left \"<\" \"<=\" \">\" \">=\";")
    LINE(    "   %left \"\\+\" \"\\-\";")
    LINE(    "   %left \"%\" \"\\*\" \"/\";")
    LINE(    "   %right UNARY_ARITHMETIC;")
    LINE(    "   %right \"!\";")
    LINE(    "   %left \"\\-\\-\";")
    LINE(    "   %left \"\\+\\+\";")
    LINE(    "   %none FUNCTION_CALL;")
    LINE(    "")
    LINE(    "   program : statements [return($0)]")
    LINE(    "           ;")
    LINE(    "   statement : expression_statement [return($0)]")
    LINE(    "             | compound_statement [return($0)]")
    LINE(    "             | forStatement [return($0)]")
    LINE(    "             | ifStatement [return($0)]")
    LINE(    "             | functionDefinition [return($0)]")
    LINE(    "             | structDefinition [return($0)]")
    LINE(    "             | \";\" [emptyStatement]")
    LINE(    "             | \"return\" expression_statement [funcReturn($1)]")
    LINE(    "             ;")
    LINE(    "   functionDefinition : \"fun\" identifier \"\\(\" variadicFuncParList \"\\)\" compound_statement [funcDef($0,$1,$3,$5)]")
    LINE(    "                      ;")
    LINE(    "   structDefinition : \"struct\" identifier \"{\" structMemberList \"}\" \";\" [structDef($0,$1,$3)]")
    LINE(    "                    ;")
    LINE(    "   structMemberList : [emptyStructMemberList]")
    LINE(    "                    | nonEmptyStructMemberList [return($0)]")
    LINE(    "                    ;")
    LINE(    "   nonEmptyStructMemberList : structMember \";\" [firstStructMemberOfList($0)]")
    LINE(    "                            | nonEmptyStructMemberList structMember \";\" [appendStructMemberToList($0,$1)]")
    LINE(    "                            ;")
    LINE(    "   structMember : identifier [structMemberCreate($0)]")
    LINE(    "                | identifier \"=\" expression [structMemberCreateWithDef($0,$2)]")
    LINE(    "                ;")
    LINE(    "   variadicFuncParList: funcParList [return($0)]")
    LINE(    "                      | variadicFuncPar [emptyVariadicParList($0)]")
    LINE(    "                      | nonEmptyFuncParList \",\" variadicFuncPar [makeVariadicParList($0,$2)]")
    LINE(    "                      ;")
    LINE(    "   variadicFuncPar: \"\\.\\.\\.\" [return(0)]")
    LINE(    "                  | \"&\" \"\\.\\.\\.\" [return(1)]")
    LINE(    "                  ;")
    LINE(    "   funcParList: [emptyFuncParList]")
    LINE(    "              | nonEmptyFuncParList [return($0)]")
    LINE(    "              ;")
    LINE(    "   nonEmptyFuncParList: funcPar [firstFuncParOfList]")
    LINE(    "                      | nonEmptyFuncParList \",\" funcPar [appendFuncParToList($0,$2)]")
    LINE(    "                       ;")
    LINE(    "   funcPar: neededFuncPar [return($0)]")
    LINE(    "          | neededFuncPar \"=\" expression [addFuncParDefault($0,$2)]")
    LINE(    "          ;")
    LINE(    "   neededFuncPar: identifier [funcParCreate(0,$0)]")
    LINE(    "                | \"&\" identifier [funcParCreate(1,$1)]")
    LINE(    "                ;")
    LINE(    "   forStatement: \"for\" \"\\(\" forInit \";\" forCheck \";\" forIncr \"\\)\" statement [forStatement($2,$4,$6,$8)]")
    LINE(    "               ;")
    LINE(    "   forInit: expression [return($0)]")
    LINE(    "          |")
    LINE(    "          ;")
    LINE(    "   forCheck: expression [return($0)]")
    LINE(    "           |")
    LINE(    "           ;")
    LINE(    "   forIncr: expression [return($0)]")
    LINE(    "          |")
    LINE(    "          ;")
    LINE(    "   ifStatement: \"if\" \"\\(\" expression \"\\)\" statement %precedence lowerThanElse [ifStatement($2,$4)]")
    LINE(    "              | \"if\" \"\\(\" expression \"\\)\" statement \"else\" statement [ifElseStatement($2,$4,$6)]")
    LINE(    "              ;")
    LINE(    "    compound_statement : \"{\" statements \"}\" [return($1)]")
    LINE(    "                       ;")
    LINE(    "    statements : [createStatements]")
    LINE(    "               | statements statement [appendStatement]")
    LINE(    "               ;")
    LINE(    "    expression_statement : expression \";\" [return($0)]")
    LINE(    "                         ;")
    LINE(    "    expression : postfix_expression %precedence \"\\+\\+\" [return($0)]")
    LINE(    "               | expression \"\\+\" expression [binarySum($0,$2)]")
    LINE(    "               | expression \"\\-\" expression [binaryDiff($0,$2)]")
    LINE(    "               | expression \"\\*\" expression [binaryProd($0,$2)]")
    LINE(    "               | expression \"/\" expression [binaryDiv($0,$2)]")
    LINE(    "               | expression \"%\" expression [binaryMod($0,$2)]")
    LINE(    "               | \"\\+\" expression %precedence UNARY_ARITHMETIC [unaryPlus($1)]")
    LINE(    "               | \"\\-\" expression %precedence UNARY_ARITHMETIC [unaryMinus($1)]")
    LINE(    "               | \"!\" expression [unaryNot($1)]")
    LINE(    "               | \"&\" expression %precedence UNARY_ARITHMETIC [unaryReference($1)]")
    LINE(    "               | expression \"<\" expression [binarySmaller($0,$2)]")
    LINE(    "               | expression \">\" expression [binaryGreater($0,$2)]")
    LINE(    "               | expression \"<=\" expression [binarySmallerEqual($0,$2)]")
    LINE(    "               | expression \">=\" expression [binaryGreaterEqual($0,$2)]")
    LINE(    "               | expression \"==\" expression [binaryCompare($0,$2)]")
    LINE(    "               | expression \"!=\" expression [binaryInequal($0,$2)]")
    LINE(    "               | expression \"&&|and\" expression [binaryAnd($0,$2)]")
    LINE(    "               | expression \"\\|\\||or\" expression [binaryOr($0,$2)]")
    LINE(    "               | assign_expression [return($0)]")
    LINE(    "               ;")
    LINE(    "    postfix_expression : primary_expression [return($0)]")
    LINE(    "                       | postfix_expression \"\\+\\+\" [unaryPostfixIncrement($0)]")
    LINE(    "                       | postfix_expression \"\\-\\-\" [unaryPostfixDecrement($0)]")
    LINE(    "                       | postfix_expression \"\\(\" funcArgsList \"\\)\" %precedence FUNCTION_CALL [funcCall($0,$2)] ")
    LINE(    "                       | postfix_expression \"\\[\" expression \"\\]\" [subscribe($0,$2)] ")
    LINE(    "                       | postfix_expression \"\\.\" identifier [memberAccess($0,$2)] ")
    LINE(    "                       ;")
    LINE(    "    funcArgsList : [emptyFuncArgList]")
    LINE(    "                 | nonemptyFuncArgsList [return($0)]")
    LINE(    "                 ;")
    LINE(    "    nonemptyFuncArgsList : arg [firstFuncArgOfList]")
    LINE(    "                         | funcArgsList \",\" arg [appendFuncArgToList($0,$2)]")
    LINE(    "                         ;")
    LINE(    "    arg : \"\\.\" identifier \"=\" expression [createFuncArg($1,$3)]")
    LINE(    "        | expression [createFuncArg(\"\",$0)]")
    LINE(    "        ;")
    LINE(    "    primary_expression : identifier [return($0)]")
    LINE(    "                       | \"[0-9]+\" [convToInt($0)]")
    LINE(    "                       | \"([0-9]+(\\.[0-9]*)?|(\\.[0-9]+))((e|E)(\\+|\\-)?[0-9]+)?\" [convToFloat($0)]")
    LINE(    "                       | \"\\\"([^\\\"\\\\]|\\\\.)*\\\"\" [convToStr($0)]")
    LINE(    "                       | \"\\(\" expression \"\\)\" %precedence BRACKETS [return($1)]")
    LINE(    "                       | \"lambda\" \"\\(\" variadicFuncParList \"\\)\" compound_statement [lambdaFuncDef($0,$2,$4)]")
    LINE(    "                       ;")
    LINE(    "    assign_expression : postfix_expression \"=\" expression %precedence \"=\" [assign($0,$2,0)]")
    LINE(    "                      | postfix_expression \"\\*=\" expression [assign($0,$2,1)]")
    LINE(    "                      | postfix_expression \"/=\" expression [assign($0,$2,2)]")
    LINE(    "                      | postfix_expression \"\\+=\" expression [assign($0,$2,3)]")
    LINE(    "                      | postfix_expression \"\\-=\" expression [assign($0,$2,4)]")
    LINE(    "                      ;")
    LINE(    "    identifier : \"[a-zA-Z_][a-zA-Z0-9_]*\" [convToId($0)]")
    LINE(    "               ;")
    LINE(    "}");
  
#undef LINE
  
#undef BARE_WHITESPACES
#undef CPP_COMMENT
#undef C_COMMENT
  
  auto grammar=pp::createGrammar(cGrammar);
  
  // Name lexical classes through their parsed semantic actions, without repeating regexes
  constexpr std::array<std::pair<std::string_view,std::string_view>,4> terminalNames{{
    {"convToId","identifier"},
    {"convToInt","integer"},
    {"convToFloat","number"},
    {"convToStr","string"}}};
  for(const auto& production : grammar.productions)
    if(production.nRhs()==1 and production.iAction)
      {
        auto& terminal=grammar.symbols[production.iRhsList.front()];
        if(terminal.type==BaseGrammarSymbol::Type::TERMINAL_SYMBOL)
          {
            const GrammarAction action(grammar.actionStrings[*production.iAction]);
            for(const auto& [name,label] : terminalNames)
              if(action.name==name)
                terminal.displayName=label;
          }
      }
  
  return grammar;
}

const ScriptLocation* scriptNodeLocation(const ASTNode& node)
{
  return std::visit([](const auto& part)->const ScriptLocation*
  {
    const ScriptLocation* location=nullptr;
    if constexpr(requires {part.location;})
      if(part.location)
        location=&*part.location;
    
    // Composite expressions reuse the first token retained by their operands
    const auto inspect=[&location](const std::shared_ptr<ASTNode>& child)
    {
      if(not location and child)
        location=scriptNodeLocation(*child);
    };
    if constexpr(requires {part.lhs;})
      inspect(part.lhs);
    else if constexpr(requires {part.arg1;})
      inspect(part.arg1);
    else if constexpr(requires {part.base;})
      inspect(part.base);
    else if constexpr(requires {part.arg;})
      inspect(part.arg);
    else if constexpr(requires {part.callee;})
      inspect(part.callee);
    else if constexpr(requires {part.subNodes;})
      for(size_t i=0;i<part.subNodes.size() and not location;i++)
        inspect(part.subNodes[i]);
    
    return location;
  },node);
}

/// Associate grammar action names with functions that construct the
/// AST The requiredActions comes from the parser; token values enter
/// through ValueNode Operator functions visit runtime alternatives
/// and reject unsupported C++ operations
pp::ParseTreeExecutor<ASTNode,ValueNode> getParseTreeExecuctor(const std::vector<std::string_view>& requiredActions,
							       const bool editor)
{
  using namespace pp::internal;
  
  using ParseTreeExecutor=
    pp::ParseTreeExecutor<ASTNode,ValueNode>;
  
  std::map<std::string,ParseTreeExecutor::ActionFun> providedActions;
  
#define ENSURE_N_SYMBOLS(NAME,N)		\
  if(subNodes.size()!=N)						\
    errorEmitter("action " NAME " expecting " #N " symbols, obtained ",subNodes.size())
  
#define PROVIDE_ACTION_WITH_N_SYMBOLS(NAME,				\
				      N,				\
				      BODY...)				\
  providedActions[NAME]=						\
    [=](std::vector<std::shared_ptr<ASTNode>>& subNodes) -> std::shared_ptr<ASTNode> \
    {									\
      ENSURE_N_SYMBOLS(NAME,N);						\
      									\
      BODY;								\
    }
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("createStatements",0,return std::make_shared<ASTNode>(ASTNodesNode{}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("emptyStatement",1,return std::make_shared<ASTNode>(ValueNode{.value{std::monostate{}}}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("appendStatement",2,
				if(ASTNodesNode* l=std::get_if<ASTNodesNode>(&*(subNodes[0]));l==nullptr)
				  errorEmitter("first argument is not a list of statement: ",
					       std::visit([](const auto& x)
					       {
						 return typeid(x).name();
					       },*(subNodes[0])));
				else
				  l->list.emplace_back(subNodes[1]);
				return subNodes[0]);
  PROVIDE_ACTION_WITH_N_SYMBOLS("convToInt",1,
    const std::string& text=unvariant<std::string>(fetch<ValueNode>(subNodes,0).value.data);
    int64_t value{};
    const auto result=std::from_chars(text.data(),text.data()+text.size(),value);
    
    const bool valid=result.ec==std::errc{} and result.ptr==text.data()+text.size();
    if(not valid and not editor)
      {
        if(const auto& location=fetch<ValueNode>(subNodes,0).location)
          location->error("integer literal outside int64_t range: ",text);
        else
          errorEmitter("integer literal outside int64_t range: ",text);
      }
    
    Value literal;
    
    if(not valid)
      literal=Value{};
    else if(value<=std::numeric_limits<int>::max())
      literal=Value{static_cast<int>(value)};
    else
      literal=Value{value};
    
    return std::make_shared<ASTNode>(ValueNode{.value=literal,.location=fetch<ValueNode>(subNodes,0).location}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("convToStr",1,
    const auto& token=fetch<ValueNode>(subNodes,0);
    return std::make_shared<ASTNode>(ValueNode{
      .value{unescapeString(std::get<std::string>(token.value.data))},.location=token.location}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("convToId",1,
    const auto& token=fetch<ValueNode>(subNodes,0);
    return std::make_shared<ASTNode>(IdNode{
      .name=std::get<std::string>(token.value.data),.description=token.description,.location=token.location}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("convToFloat",1,
    const auto& token=fetch<ValueNode>(subNodes,0);
    return std::make_shared<ASTNode>(ValueNode{
      .value{strtod(std::get<std::string>(token.value.data).c_str(),nullptr)},.location=token.location}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("return",1,return subNodes[0]);
  PROVIDE_ACTION_WITH_N_SYMBOLS("createFuncArg",2,
				std::string name{};
				std::optional<ScriptLocation> location;
				if(IdNode* id=std::get_if<IdNode>(&*subNodes[0]))
				  {
				    name=id->name;
				    location=id->location;
				  }
				return std::make_shared<ASTNode>(FuncArgNode{.name=name,.expr=subNodes[1],.location=location}));
  
#define DEFINE_UNOP(OP,NAME)						\
  PROVIDE_ACTION_WITH_N_SYMBOLS("unary" #NAME,1,return std::make_shared<ASTNode>(UnOpNode{.arg=subNodes[0],.kind=UnOpNode::Kind::NAME, \
	  .op=[](const Value& arg)					\
	  {								\
	    return std::visit([](const auto& arg) -> Value		\
	    {								\
	      if constexpr(requires {OP arg;})				\
		return {OP arg};					\
	      else							\
		{							\
		  errorEmitter("Cannot " #NAME" the type: ",typeid(arg).name()); \
		  return {std::monostate{}};				\
		}},							\
			      arg.data);				\
	  }}))
  
  DEFINE_UNOP(!,Not);
  DEFINE_UNOP(+,Plus);
  DEFINE_UNOP(-,Minus);
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("unaryReference",1,return std::make_shared<ASTNode>(UnaryRefNode{.arg=subNodes[0]}));
  
#undef DEFINE_UNOP
  
#define DEFINE_BINOP(OP,NAME)						\
  PROVIDE_ACTION_WITH_N_SYMBOLS("binary" #NAME,2,return std::make_shared<ASTNode>(BinOpNode{.arg1=subNodes[0], \
	  .arg2=subNodes[1],						\
	  .op=[](const Value& arg1,					\
		 const Value& arg2)					\
	  {								\
	    return std::visit([](const auto& arg1,			\
				 const auto& arg2) -> Value		\
	    {								\
	      if constexpr(requires {arg1 OP arg2;})			\
		return {arg1 OP arg2};					\
	      else							\
		{							\
		  errorEmitter("Cannot " #NAME" the types: ",typeid(arg1).name()," and ",typeid(arg2).name()); \
		  return {std::monostate{}};				\
		}},							\
			      arg1.data,arg2.data);			\
	  }}))
  
#define DEFINE_REGISTERED_BINOP(NAME,OPERATOR) \
  PROVIDE_ACTION_WITH_N_SYMBOLS("binary" #NAME,2, \
    return std::make_shared<ASTNode>(BinOpNode{ \
      .arg1=subNodes[0],.arg2=subNodes[1],.operatorType=binaryOperatorType<OPERATOR>(),.op=applyBinaryOperator<OPERATOR>}))
  
  DEFINE_REGISTERED_BINOP(Sum,SumOperator);
  DEFINE_REGISTERED_BINOP(Diff,DifferenceOperator);
  DEFINE_REGISTERED_BINOP(Prod,ProductOperator);
  DEFINE_REGISTERED_BINOP(Div,DivisionOperator);
  DEFINE_REGISTERED_BINOP(Mod,RemainderOperator);
  
#undef DEFINE_REGISTERED_BINOP
  
  DEFINE_BINOP(<,Smaller);
  DEFINE_BINOP(>,Greater);
  DEFINE_BINOP(<=,SmallerEqual);
  DEFINE_BINOP(>=,GreaterEqual);
  DEFINE_BINOP(==,Compare);
  DEFINE_BINOP(!=,Inequal);
  
#undef DEFINE_BINOP
  
#define DEFINE_LOGICAL_BINOP(NAME,MODE) \
  PROVIDE_ACTION_WITH_N_SYMBOLS("binary" #NAME,2, \
    return std::make_shared<ASTNode>(BinOpNode{ \
      .arg1=subNodes[0],.arg2=subNodes[1], \
      .evaluationMode=BinOpNode::EvaluationMode::MODE}))
  
  DEFINE_LOGICAL_BINOP(Or,LOGICAL_OR);
  DEFINE_LOGICAL_BINOP(And,LOGICAL_AND);
  
#undef DEFINE_LOGICAL_BINOP
  
#define DEFINE_PREPOSTFIX_OP(NAME,PRE,POST)					\
  PROVIDE_ACTION_WITH_N_SYMBOLS("unary" #NAME,1,return std::make_shared<ASTNode>(PrePostfixOpNode{.arg=subNodes[0], \
	  .op=[](Value& arg)						\
	  {								\
	    return std::visit([](auto& arg) -> Value			\
	    {								\
	      if constexpr(requires {PRE arg POST;})			\
		return {PRE arg POST};					\
	      else							\
		{							\
		  errorEmitter("Cannot " #NAME" the type: %s",typeid(arg).name()); \
		  return {std::monostate{}};				\
		}},							\
			      arg.data);				\
	  }}))
  
  DEFINE_PREPOSTFIX_OP(PostfixIncrement,,++);
  //  DEFINE_PREPOSTFIX_OP(PrefixIncrement,++,);
  DEFINE_PREPOSTFIX_OP(PostfixDecrement,,--);
  //DEFINE_PREPOSTFIX_OP(PrefixDecrement,--,);
  
#undef DEFINE_PREPOSTFIX_OP
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("structMemberCreate",1,
    const auto& id=fetch<IdNode>(subNodes,0);
    return std::make_shared<ASTNode>(StructMemberNode{.name=id.name,.location=id.location,.description=id.description}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("structMemberCreateWithDef",2,
    const auto& id=fetch<IdNode>(subNodes,0);
    return std::make_shared<ASTNode>(StructMemberNode{.name=id.name,.def=subNodes[1],.location=id.location,.description=id.description}));
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("funcParCreate",2,
    const auto& id=fetch<IdNode>(subNodes,1);
    return std::make_shared<ASTNode>(FuncParNode{
      .name=id.name,.isRef=(bool)unvariant<int>(fetch<ValueNode>(subNodes,0).value.data),.location=id.location}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("addFuncParDefault",2,
    const FuncParNode& in=fetch<FuncParNode>(subNodes,0);
    return std::make_shared<ASTNode>(FuncParNode{.name=in.name,.isRef=in.isRef,.def=subNodes[1],.location=in.location}));
#define PROVIDE_FUNC_LIST_ACTIONS(NAME)					\
  PROVIDE_ACTION_WITH_N_SYMBOLS("empty" #NAME "List",0,return std::make_shared<ASTNode>(NAME ## ListNode{})); \
  PROVIDE_ACTION_WITH_N_SYMBOLS("first" #NAME "OfList",1,return makeFirstOfList<NAME ## Node,NAME ## ListNode>(subNodes[0])); \
  PROVIDE_ACTION_WITH_N_SYMBOLS("append" #NAME "ToList",2,return appendToList<NAME ## Node,NAME ## ListNode>(subNodes[0],subNodes[1])); \
  
  PROVIDE_FUNC_LIST_ACTIONS(FuncPar);
  PROVIDE_FUNC_LIST_ACTIONS(FuncArg);
  PROVIDE_FUNC_LIST_ACTIONS(StructMember);
  
#undef PROVIDE_FUNC_LIST_ACTIONS
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("emptyVariadicParList",1,return std::make_shared<ASTNode>(FuncParListNode::getEmptyVariadic(unvariant<int>(fetch<ValueNode>(subNodes,0).value.data))));
  PROVIDE_ACTION_WITH_N_SYMBOLS("makeVariadicParList",2,
				fetch<FuncParListNode>(subNodes,0).makeVariadic(unvariant<int>(fetch<ValueNode>(subNodes,1).value.data));
				return subNodes[0]);
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("funcDef",4,
    return std::make_shared<ASTNode>(FuncDefNode{
      .name=fetch<IdNode>(subNodes,1).name,
      .fun=std::make_shared<FuncNode>(fetch<FuncParListNode>(subNodes,2),subNodes[3],
                                     fetch<ValueNode>(subNodes,0).description,
                                     fetch<ValueNode>(subNodes,0).location),
      .location=fetch<IdNode>(subNodes,1).location}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("structDef",3,
    return std::make_shared<ASTNode>(StructDefNode{
      .name=fetch<IdNode>(subNodes,1).name,
      .members=fetch<StructMemberListNode>(subNodes,2),
      .description=fetch<ValueNode>(subNodes,0).description,
      .location=fetch<ValueNode>(subNodes,0).location}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("lambdaFuncDef",3,
    return std::make_shared<ASTNode>(FuncNode{
      fetch<FuncParListNode>(subNodes,1),subNodes[2],fetch<ValueNode>(subNodes,0).description,
      fetch<ValueNode>(subNodes,0).location}));
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("funcCall",2,
    if(const auto* id=std::get_if<IdNode>(subNodes[0].get()))
      return std::make_shared<ASTNode>(FuncCallNode{
        .fun=id->name,.args=fetch<FuncArgListNode>(subNodes,1),.location=id->location});
    else
      return std::make_shared<ASTNode>(FuncCallNode{
        .args=fetch<FuncArgListNode>(subNodes,1),.callee=subNodes[0]}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("funcReturn",1,return std::make_shared<ASTNode>(ReturnNode{.arg=subNodes[0]}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("subscribe",2,return std::make_shared<ASTNode>(SubscribeNode{.base=subNodes[0],.subscr=subNodes[1]}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("memberAccess",2,return std::make_shared<ASTNode>(MemberAccessNode{
    .base=subNodes[0],.member=fetch<IdNode>(subNodes,1).name,.memberLocation=fetch<IdNode>(subNodes,1).location}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("forStatement",4,return std::make_shared<ASTNode>(ForNode{.subNodes{subNodes}}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("ifStatement",2,return std::make_shared<ASTNode>(IfNode{.subNodes{subNodes}}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("ifElseStatement",3,return std::make_shared<ASTNode>(IfNode{.subNodes{subNodes}}));
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("assign",3,
    const AssignMode mode=static_cast<AssignMode>(unvariant<int>(fetch<ValueNode>(subNodes,2).value.data));
    
    // Documentation before a named lambda binding belongs to the callable as well
    if(mode==AssignMode::PLAIN)
      if(const auto* id=std::get_if<IdNode>(subNodes[0].get());id and not id->description.empty())
        if(auto* function=std::get_if<FuncNode>(subNodes[1].get());function and function->description.empty())
          function->description=id->description;
    
    return std::make_shared<ASTNode>(AssignNode{.lhs=subNodes[0],.rhs=subNodes[1],.mode=mode}));
  
#undef PROVIDE_ACTION_WITH_N_SYMBOLS
#undef ENSURE_N_SYMBOLS
  
  return ParseTreeExecutor(providedActions,requiredActions);
}

/// Register built-in functions and constants in the evaluator root scope
void registerEvaluatorBuiltins(Evaluator& evaluator)
{
  Environment& env=*evaluator.curEnv;
  
  env["M_PI"]=M_PI;
  env.describe("M_PI","Value of pi matching C definition");
  
  const auto registerListing=[&env](const char* name,
                                   const EntitySelection selection,
                                   const char* description)
  {
    env[name]=HostFunction{
      [selection,name](std::vector<std::shared_ptr<Value>>& args,
                       Evaluator& caller)->Value
      {
        const int verbose=args.empty() or not args[0]?1:
          hostArgumentTypeConcretizer<int>(args,0,name);
	
        if(args.size()>1 or (verbose!=0 and verbose!=1))
          errorEmitter(name," accepts only an optional verbose argument, 0 or 1");
	
        const std::string text=formatEntityListing(*caller.curEnv,selection,verbose!=0);
        if(not text.empty())
          MASTER_PRINTF("%s",text.c_str());
	
        return {};
      },{"verbose"},false,description,{{"[verbose = 1]"}}};
  };
  registerListing("listFunctions",EntitySelection::Functions,"List visible callables (default: verbose=1)");
  registerListing("listVariables",EntitySelection::Variables,"List visible non-callable bindings (default: verbose=1)");
  registerListing("listEntities",EntitySelection::All,"List all visible bindings (default: verbose=1)");
  
  env["help"]=HostFunction{
    [](std::vector<std::shared_ptr<Value>>& args,
       Evaluator& caller)->Value
    {
      if(args.size()>1)
        errorEmitter("help accepts only an optional keyword string");
      
      const std::string keyword=args.empty() or not args[0]?"":
        hostArgumentTypeConcretizer<std::string>(args,0,"help");
      const std::string text=formatHelp(*caller.curEnv,keyword);
      MASTER_PRINTF("%s",text.c_str());
      
      return {};
    },{"keyword"},false,
    "Explain the script language, or a named entity visible in the calling scope (default: keyword=\"\")",
    {{"[keyword = \"\"]"}}};
  
  registerHostFunction(evaluator,"describe",{"name","description"},
                       [](Evaluator& caller,
                          const std::string& name,
                          const std::string& description)
  {
    caller.curEnv->describe(name,description);
  },"Attach a description to a visible binding");
  
  env["print"]=
    HostFunction{
    [](std::vector<std::shared_ptr<Value>>& args,
	 Evaluator&)->Value
    {
	std::ostringstream output;
	for(std::shared_ptr<Value>& arg : args)
	  std::visit([&output](const auto& v)
	  {
	    if constexpr(requires {output<<v;})
	      output<<v;
	    else
	      output<<"unprintable type: "<<typeid(v).name()<<"\n";
	  },arg->data);
	
	const std::string text=output.str();
	if(not text.empty())
	  MASTER_PRINTF("%s",text.c_str());
	
	return {};
    },{},false,"Print arguments on the ranks selected by the library printing policy; no newline is added",{{"..."}}};
  
  env["makeVec"]=
    HostFunction{
    [](std::vector<std::shared_ptr<Value>>& args,
	 Evaluator&)->Value
    {
	return {ValuesList{.data{args}}};
    },{},false,"Create a list from the supplied arguments",{{"..."}}};
  
  env["forEachEl"]=
    HostFunction{
    [](std::vector<std::shared_ptr<Value>>& args,
	 Evaluator& ev)->Value
    {
	if(const size_t n=args.size();n!=2)
	  errorEmitter("trying to call forEachEl with ",n," elements in place of 2");
	
	std::vector<std::shared_ptr<Value>> v(1);
	if(ValuesList* list=std::get_if<ValuesList>(&args[0]->data))
	  for(std::shared_ptr<Value> val : list->data)
	    {
	      v[0]=val;
	      
	      ev.callFunction(*args[1],v);
	    }
	else
	  errorEmitter("trying to call \"forEachEl\" on a non-vector type ",variantInnerTypeName(args[0]->data));
	
	return {};
    },{},false,"Call a function on each element of a list",{{"values","function"}}};
  
  env["seq"]=
    HostFunction{
    [](std::vector<std::shared_ptr<Value>>& args,
	 Evaluator&)->Value
    {
	auto p=
	  [&args](const size_t& i)->int*
	  {
	    if(i>=args.size())
	      return nullptr;
	    else
	      {
		int* r=std::get_if<int>(&args[i]->data);
		if(not r)
		  errorEmitter("arg ",i," is of type ",variantInnerTypeName(args[i]->data)," instead of ",typeid(int).name()," (int)");
		
		return r;
	      }
	  };
	
	if(p(3))
	  errorEmitter("too many argument to init the vec of int");
	
	int min{0};
	int delta{1};
	int max{};
	if(int* p0=p(0))
	  if(int* p1=p(1))
	    if(int* p2=p(2))
	      {
		max=*p2;
		delta=*p1;
		min=*p0;
	      }
	    else
	      {
		max=*p1;
		min=*p0;
	      }
	  else
	    max=*p0;
	else
	  errorEmitter("too few arguments to init a vec of int, passed ",args.size());
	
	if(delta==0)
	  errorEmitter("asked to create a sequence with zero step");
	
	ValuesList v;
	for(int i=min;i<max;i+=delta)
	  v.data.push_back(std::make_shared<Value>(i));
	
	return {v};
    },{},false,"Create integers below max: seq(max), seq(min,max), or seq(min,step,max); omitted min=0 and step=1",{{"max"},{"min","max"},{"min","step","max"}}};
  
  const auto requireArgCount=[](const std::vector<std::shared_ptr<Value>>& args,
                                const size_t expected,
                                const char* name)
  {
    if(args.size()!=expected)
      errorEmitter("trying to call function ",name," which expects ",expected,
                   " args with ",args.size());
  };
  
  const auto numericArg=[](const Value& value,
                           const char* name)->double
  {
    if(const int* integer=std::get_if<int>(&value.data))
      return *integer;
    else if(const int64_t* integer=std::get_if<int64_t>(&value.data))
      return *integer;
    else if(const double* real=std::get_if<double>(&value.data))
      return *real;
    else
      errorEmitter("Trying to call ",name," with a nonnumeric argument");
  };
  
  env["rand"]=
    HostFunction{[requireArgCount](std::vector<std::shared_ptr<Value>>& args,
                                   Evaluator&)->Value
    {
      requireArgCount(args,0,"rand");
      
      return {std::rand()};
    },{},false,"Draw an integer from the C library RNG, independently of the simulation RNG",{{}}};
  
  env["srand"]=
    HostFunction{[requireArgCount,numericArg](std::vector<std::shared_ptr<Value>>& args,
                                              Evaluator&)->Value
    {
      requireArgCount(args,1,"srand");
      std::srand(static_cast<unsigned int>(numericArg(*args[0],"srand")));
      
      return {std::monostate{}};
    },{},false,"Seed the C library RNG without changing the simulation RNG",{{"seed"}}};
  
  const auto registerUnaryMath=[&env,requireArgCount,numericArg](const char* name,
                                                                 double (*function)(double),
                                                                 const char* description)
  {
    env[name]=
      HostFunction{[name,function,requireArgCount,numericArg](std::vector<std::shared_ptr<Value>>& args,
                                                               Evaluator&)->Value
      {
        requireArgCount(args,1,name);
	
        return {function(numericArg(*args[0],name))};
      },{},false,description,{{"value"}}};
  };
  
  registerUnaryMath("sqrt",static_cast<double (*)(double)>(std::sqrt),"Return the square root");
  registerUnaryMath("exp",static_cast<double (*)(double)>(std::exp),"Return the exponential");
  registerUnaryMath("sin",static_cast<double (*)(double)>(std::sin),"Return the sine of an angle in radians");
  registerUnaryMath("cos",static_cast<double (*)(double)>(std::cos),"Return the cosine of an angle in radians");
  registerUnaryMath("tan",static_cast<double (*)(double)>(std::tan),"Return the tangent of an angle in radians");
  
  env["pow"]=
    HostFunction{[requireArgCount,numericArg](std::vector<std::shared_ptr<Value>>& args,
                                              Evaluator&)->Value
    {
      requireArgCount(args,2,"pow");
      
      return {std::pow(numericArg(*args[0],"pow"),numericArg(*args[1],"pow"))};
    },{},false,"Raise a numeric base to a numeric exponent",{{"base","exponent"}}};
  
  env["to_string"]=
    HostFunction{[requireArgCount](std::vector<std::shared_ptr<Value>>& args,
                                   Evaluator&)->Value
    {
      requireArgCount(args,1,"to_string");
      
      if(const int* integer=std::get_if<int>(&args[0]->data))
        return {std::to_string(*integer)};
      else if(const int64_t* integer=std::get_if<int64_t>(&args[0]->data))
        return {std::to_string(*integer)};
      else if(const double* real=std::get_if<double>(&args[0]->data))
        return {std::to_string(*real)};
      else
        errorEmitter("Trying to call to_string with a nonnumeric argument");
    },{},false,"Convert a numeric value to a string",{{"value"}}};
}
