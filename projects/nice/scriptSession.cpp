#include "scriptSession.hpp"
#include "scriptAnalysis.hpp"

#include <charconv>
#include <limits>
#include <sstream>

namespace
{
  /// Parse a complete unsigned protocol field without accepting signs or trailing text
  bool readSize(std::istream& input,
                size_t& value)
  {
    std::string token;
    input>>token;
    const auto parsed=std::from_chars(token.data(),token.data()+token.size(),value);
    
    return input and parsed.ec==std::errc{} and parsed.ptr==token.data()+token.size();
  }
  
  /// Read exactly the announced bytes; source contents never act as protocol delimiters
  bool readBytes(std::istream& input,
                 const size_t size,
                 std::string& value)
  {
    const bool valid=size<=value.max_size() and size<=static_cast<size_t>(std::numeric_limits<std::streamsize>::max());
    if(valid)
      {
        value.resize(size);
        input.read(value.data(),static_cast<std::streamsize>(size));
      }
    
    return valid and input;
  }
}

int runScriptSession(const Environment& environment,
                     std::istream& input,
                     std::ostream& output,
                     std::ostream& diagnostic)
{
  ScriptAnalysisSession session;
  std::string header;
  bool valid=true;
  while(valid and std::getline(input,header))
    {
      std::istringstream fields(header);
      std::string version,operation;
      size_t id=0,offset=0,nPieces=0;
      fields>>version;
      valid=version=="NICE/1" and readSize(fields,id);
      fields>>operation;
      valid=valid and readSize(fields,offset) and readSize(fields,nPieces) and nPieces>0;
      fields>>std::ws;
      valid=valid and fields.eof();
      ScriptSource source;
      size_t begin=0,lastSize=0;
      for(size_t i=0;valid and i<nPieces;i++)
        {
          valid=static_cast<bool>(std::getline(input,header));
          std::istringstream lengths(header);
          size_t nameSize=0,textSize=0;
          valid=valid and readSize(lengths,nameSize) and readSize(lengths,textSize);
          lengths>>std::ws;
          valid=valid and lengths.eof();
          std::string name,text;
          valid=valid and readBytes(input,nameSize,name) and readBytes(input,textSize,text);
          if(valid)
            {
              begin=source.text.size();
              lastSize=text.size();
              source.append(text,name);
            }
        }
      if(valid)
        {
          const bool accepted=(operation=="check" or operation=="complete") and offset<=lastSize;
          output<<'('<<id<<' '<<(accepted?0:1)<<' ';
          std::string payload="nil";
          if(accepted)
            {
              session.update(source);
              if(operation=="check")
                payload=formatScriptDiagnosticsForEmacs(*session.parsed);
              else
                {
                  auto completion=session.complete(environment,begin+offset);
                  if(completion.applicable)
                    {
                      completion.begin-=begin;
                      completion.end-=begin;
                    }
                  payload=formatScriptCompletionForEmacs(completion);
                }
            }
          output<<session.revision<<' '<<session.nParses<<' '<<payload<<' ';
          writeScriptLispString(output,accepted?"":"Invalid editor operation or cursor offset");
          output<<")\n"<<std::flush;
          valid=static_cast<bool>(output);
        }
    }
  if(not valid or input.bad())
    diagnostic<<"NICE editor session: invalid or truncated request/response\n";
  
  return valid and not input.bad()?0:1;
}
