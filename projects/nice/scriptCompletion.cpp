#include "scriptCompletion.hpp"

#include <sstream>

namespace
{
  /// Match the ASCII identifier alphabet used by the grammar
  bool identifierByte(const char c)
  {
    return (c>='a' and c<='z') or (c>='A' and c<='Z') or
      (c>='0' and c<='9') or c=='_';
  }
  
  /// Escape text as a Lisp string, including control bytes
  void writeLispString(std::ostream& output,
                       const std::string& value)
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
}

ScriptCompletion completeScriptAt(const Environment& environment,
                                  const std::string_view text,
                                  const size_t byteOffset)
{
  ScriptCompletion result;
  
  if(byteOffset<=text.size())
    {
      result.begin=byteOffset;
      result.end=byteOffset;
      while(result.begin>0 and identifierByte(text[result.begin-1]))
        result.begin--;
      while(result.end<text.size() and identifierByte(text[result.end]))
        result.end++;
      
      // Scan only the prefix, accepting unfinished calls and rejecting comments and strings
      enum class State{Code,String,LineComment,BlockComment};
      State state=State::Code;
      bool escaped=false;
      std::vector<std::string_view> tokens;
      std::vector<size_t> delimiters;
      for(size_t i=0;i<result.begin;i++)
        {
          const char c=text[i];
          const char next=i+1<result.begin?text[i+1]:'\0';
          
          if(state==State::String)
            {
              if(c=='"' and not escaped)
                state=State::Code;
              
              escaped=c=='\\' and not escaped;
            }
          else if(state==State::LineComment)
            {
              if(c=='\n' or c=='\r')
                state=State::Code;
            }
          else if(state==State::BlockComment)
            {
              if(c=='*' and next=='/')
                {
                  state=State::Code;
                  i++;
                }
            }
          else if(c=='/' and (next=='/' or next=='*'))
            {
              state=next=='/'?State::LineComment:State::BlockComment;
              i++;
            }
          else if(c=='"')
            {
              state=State::String;
              tokens.push_back(text.substr(i,1));
            }
          else if(identifierByte(c))
            {
              const size_t begin=i;
              while(i+1<result.begin and identifierByte(text[i+1]))
                i++;
              
              tokens.push_back(text.substr(begin,i+1-begin));
            }
          else if(c!=' ' and c!='\t' and c!='\r' and c!='\n' and c!='\f' and c!='\v')
            {
              if(c=='(' or c=='[' or c=='{')
                delimiters.push_back(tokens.size());
              else if((c==')' or c==']' or c=='}') and not delimiters.empty())
                delimiters.pop_back();
              
              tokens.push_back(text.substr(i,1));
            }
        }
      
      const auto prefix=text.substr(result.begin,byteOffset-result.begin);
      const bool numeric=not prefix.empty() and prefix.front()>='0' and prefix.front()<='9';
      if(state==State::Code and not numeric)
        {
          result.applicable=true;
          const auto add=[&](const std::string& name,
                            const EntityInfo& info,
                            const std::string& documentation)
          {
            if(name.starts_with(prefix))
              result.candidates.push_back({name,info,documentation});
          };
          const auto names=environment.entityNames();
          const auto infoOf=[&](const std::string_view name)->std::optional<EntityInfo>
          {
            const std::string key(name);
            std::optional<EntityInfo> info;
            if(std::binary_search(names.begin(),names.end(),key))
              info=environment.entityInfo(key);
            
            return info;
          };
          
          if(not tokens.empty() and tokens.back()==".")
            {
              // Complete f(.parameter or f(value,.parameter from its registered signature
              if(tokens.size()>1 and (tokens[tokens.size()-2]=="(" or tokens[tokens.size()-2]==","))
                {
                  if(not delimiters.empty())
                    {
                      const size_t open=delimiters.back();
                      if(open>0 and tokens[open]=="(" and (open<2 or tokens[open-2]!="."))
                        if(const auto info=infoOf(tokens[open-1]);info and info->callable)
                          for(const auto& form : info->argumentForms)
                            for(auto argument : form)
                              {
                                const size_t first=argument.find_first_not_of("[& ");
                                argument=first==std::string::npos?"":argument.substr(first);
                                argument=argument.substr(0,argument.find_first_of(" =]"));
                                if(not argument.empty() and argument!="..." and argument.front()!='<')
                                  add(argument,EntityInfo{"parameter"},formatHelp(environment,std::string(tokens[open-1])));
                              }
                    }
                }
              else if(tokens.size()>1 and (tokens.size()<3 or tokens[tokens.size()-3]!="."))
                {
                  // Inspect storage metadata only; never evaluate member access or getters
                  const std::string objectName(tokens[tokens.size()-2]);
                  const Environment* scope=&environment;
                  while(scope and not scope->varTable.contains(objectName))
                    scope=scope->parent.get();
                  
                  if(scope and scope->varTable.at(objectName))
                    {
                      const auto& value=scope->varTable.at(objectName)->data;
                      if(const auto* object=std::get_if<HObject>(&value);object and object->info)
                        {
                          for(const auto& [name,accessor] : object->info->memberAccessorsProvider)
                            add(name,EntityInfo{"host member"},"");
                        }
                      else if(const auto* object=std::get_if<Struct>(&value);object and object->structDefNode)
                        for(const auto& member : object->structDefNode->members.list)
                          add(member.name,EntityInfo{.kind="script member",.location=member.location},"");
                    }
                }
            }
          else
            {
              for(const auto& name : names)
                add(name,environment.entityInfo(name),formatHelp(environment,name));
              
              for(const std::string keyword : {"fun","lambda","struct","if","else","for","return","and","or"})
                if(not std::binary_search(names.begin(),names.end(),keyword))
                  add(keyword,EntityInfo{"keyword"},"");
            }
          
          std::sort(result.candidates.begin(),result.candidates.end(),
            [](const auto& lhs,
               const auto& rhs)
            {
              return lhs.name<rhs.name;
            });
          result.candidates.erase(std::unique(result.candidates.begin(),result.candidates.end(),
            [](const auto& lhs,
               const auto& rhs)
            {
              return lhs.name==rhs.name;
            }),result.candidates.end());
        }
    }
  
  return result;
}

std::string formatScriptCompletionForEmacs(const ScriptCompletion& completion)
{
  std::ostringstream output;
  if(not completion.applicable)
    output<<"nil";
  else
    {
      output<<'('<<completion.begin<<' '<<completion.end<<" (";
      for(const auto& candidate : completion.candidates)
        {
          output<<'(';
          writeLispString(output,candidate.name);
          output<<' ';
          writeLispString(output,candidate.info.kind);
          output<<' ';
          writeLispString(output,candidate.documentation);
          output<<' ';
          writeLispString(output,candidate.info.location?candidate.info.location->source:"");
          const SourceCursor cursor=candidate.info.location?candidate.info.location->cursor:SourceCursor{};
          output<<' '<<cursor.line<<' '<<cursor.column<<' '<<cursor.offset<<')';
        }
      output<<"))";
    }
  
  return output.str();
}
