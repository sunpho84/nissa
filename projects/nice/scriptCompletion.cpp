#include "scriptCompletion.hpp"

#include <sstream>
#include <map>
#include <span>
#include <charconv>

namespace
{
  /// Match the ASCII identifier alphabet used by the grammar
  bool identifierByte(const char c)
  {
    return (c>='a' and c<='z') or (c>='A' and c<='Z') or
      (c>='0' and c<='9') or c=='_';
  }
  
  struct CompletionFunction;
  
  /// Metadata for a visible value, without constructing or reading its payload
  struct CompletionBinding
  {
    /// Classification and declaration position offered to clients
    EntityInfo info{"script variable (unknown type)"};
    
    /// Help retained for existing registered entities and callable aliases
    std::string documentation;
    
    /// Members of this value, when its type is known
    std::map<std::string,EntityInfo> members;
    
    /// Concrete scalar or host payload identity, when known
    std::optional<std::type_index> cppType{};
    
    /// Host metadata used to select binary overloads
    std::shared_ptr<const HObjectInfo> hostType{};
    
    /// Result of calling this binding, absent when its result is dynamic or unknown
    std::shared_ptr<CompletionBinding> callResult;
    
    /// Unevaluated body of a script factory, borrowed only during this query
    std::shared_ptr<CompletionFunction> scriptFunction;
  };
  
  /// Restricted script body whose local return expressions can be typed without execution
  struct CompletionFunction
  {
    /// Complete body tokens, excluding its enclosing braces
    std::span<const std::string_view> body;
    
    /// Fixed host callables visible at definition, checked again at the call site
    std::map<std::string,std::shared_ptr<CompletionBinding>> hostResults;
  };
  
  /// Infer agreeing reachable returns; parameters and captured values remain unknown
  CompletionBinding scriptReturn(const CompletionFunction& function,
                                 const std::map<std::string,CompletionBinding>& bindings);
  
  /// Collect names from shared host metadata without binding or calling accessors
  std::map<std::string,EntityInfo> hostMembers(const HObjectInfo& type)
  {
    std::map<std::string,EntityInfo> members;
    for(const auto& [name,provider] : type.memberAccessorsProvider)
      members.emplace(name,EntityInfo{"host member"});
    
    return members;
  }
  
  /// Collect declared script members without evaluating their defaults
  std::map<std::string,EntityInfo> scriptMembers(const StructDefNode& type)
  {
    std::map<std::string,EntityInfo> members;
    for(const auto& member : type.members.list)
      members.emplace(member.name,EntityInfo{.kind="script member",.location=member.location});
    
    return members;
  }
  
  /// Recognize a complete identifier rather than a numeric token or punctuation
  bool identifier(const std::string_view token)
  {
    return not token.empty() and not (token.front()>='0' and token.front()<='9') and
      std::all_of(token.begin(),token.end(),identifierByte);
  }
  
  /// Locate a token in the supplied editor text
  ScriptLocation completionLocation(const ScriptSource& source,
                                    const std::string_view token)
  {
    const auto location=source.locate(SourceCursor{.offset=static_cast<size_t>(token.data()-source.text.data())});
    
    return {std::string(location.source),location.cursor};
  }
  
  /// Reuse the language's immediately preceding documentation-comment rules
  std::string completionDescription(const ScriptSource& source,
                                   const std::span<const std::string_view> tokens,
                                   const size_t index)
  {
    const size_t begin=index==0?0:tokens[index-1].data()+tokens[index-1].size()-source.text.data();
    const size_t end=tokens[index].data()-source.text.data();
    const auto cursor=source.locate(SourceCursor{.offset=begin}).cursor;
    
    return precedingDocumentation(pp::internal::ShiftNode{
      .txt=tokens[index],.leadingTrivia=std::string_view(source.text).substr(begin,end-begin),.triviaCursor=cursor});
  }
  
  /// Format a statically recovered signature and its declaration documentation
  std::string completionDocumentation(const std::string& name,
                                      const EntityInfo& info)
  {
    std::ostringstream output;
    for(const auto& form : info.argumentForms)
      {
        output<<name<<'(';
        for(size_t i=0;i<form.size();i++)
          output<<(i==0?"":", ")<<form[i];
        output<<")\n";
      }
    
    if(not info.description.empty())
      output<<info.description<<'\n';
    if(info.location)
      output<<info.location->describe()<<'\n';
    
    return output.str();
  }
  
  /// One syntactic parameter; defaults describe optionality, never an assumed value type
  struct CompletionParameter
  {
    /// Name visible in the function body, including the reserved variadic binding
    std::string name;
    
    /// Signature text retaining references, optionality, and the unevaluated default
    std::string form;
    
    /// Parameter declaration position
    ScriptLocation location;
  };
  
  /// Recover comma-separated parameters, respecting nested calls, lists, and string tokens
  std::vector<CompletionParameter> completionParameters(const std::span<const std::string_view> tokens,
                                                        const ScriptSource& source)
  {
    std::vector<CompletionParameter> result;
    size_t begin=0;
    size_t depth=0;
    for(size_t i=0;i<=tokens.size();i++)
      {
        if(i==tokens.size() or (tokens[i]=="," and depth==0))
          {
            const auto parameter=tokens.subspan(begin,i-begin);
            const bool reference=not parameter.empty() and parameter.front()=="&";
            const size_t nameIndex=reference?1:0;
            if(nameIndex<parameter.size())
              {
                const bool variadic=parameter.size()==nameIndex+3 and
                  parameter[nameIndex]=="." and parameter[nameIndex+1]=="." and parameter[nameIndex+2]==".";
                if(variadic or identifier(parameter[nameIndex]))
                  {
                    const std::string name=variadic?FuncParNode::VARIADIC_PAR_NAME:std::string(parameter[nameIndex]);
                    std::string form=(reference?"&":"")+(variadic?std::string("..."):name);
                    if(parameter.size()>nameIndex+2 and parameter[nameIndex+1]=="=")
                      {
                        const auto first=parameter[nameIndex+2];
                        const auto last=parameter.back();
                        const std::string_view expression(first.data(),last.data()+last.size()-first.data());
                        form="["+form+" = "+std::string(expression)+"]";
                      }
                    
                    result.push_back({name,std::move(form),completionLocation(source,parameter[nameIndex])});
                  }
              }
            
            begin=i+1;
          }
        else if(tokens[i]=="(" or tokens[i]=="[")
          depth++;
        else if((tokens[i]==")" or tokens[i]=="]") and depth>0)
          depth--;
      }
    
    return result;
  }
  
  /// Static scope reconstructed from the token prefix, separate from runtime storage
  struct CompletionScope
  {
    /// Visible names, copied on entry so deferred bodies cannot change their parents
    std::map<std::string,CompletionBinding> bindings;
    
    /// Names assigned in this scope, including writes to enclosing bindings
    std::set<std::string> written;
    
    /// First token of the current statement
    size_t begin{0};
    
    /// Enclosing expression depth restored after this block
    size_t parentDepth{0};
    
    /// Plain blocks propagate assignments; conditional blocks discard inferred result types
    enum class Kind{Plain,Conditional,Deferred};
    
    /// Whether this block executes directly, conditionally, or as a deferred function body
    Kind kind{Kind::Plain};
    
    /// Named definition or direct lambda assignment published when its body closes
    std::string definitionName;
    
    /// Callable metadata under construction for this deferred body
    std::optional<CompletionBinding> definition;
    
    /// Whether semicolon-terminated statements declare struct members
    bool structDefinition{false};
    
    /// A direct lambda initializer must end at its body to remain a callable alias
    bool lambdaDefinition{false};
  };
  
  /// Build a value description from static metadata, without constructing its payload
  CompletionBinding resultBinding(const ValueTypeInfo& type)
  {
    CompletionBinding result;
    result.cppType=type.cppType;
    result.hostType=type.hostType.lock();
    if(result.hostType)
      {
        result.info.kind="host object '"+result.hostType->name+"'";
        result.members=hostMembers(*result.hostType);
      }
    else if(type.cppType)
      {
        if(*type.cppType==typeid(int))
          result.info.kind="variable (int)";
        else if(*type.cppType==typeid(int64_t))
          result.info.kind="variable (int64)";
        else if(*type.cppType==typeid(double))
          result.info.kind="variable (double)";
        else if(*type.cppType==typeid(std::string))
          result.info.kind="variable (string)";
      }
    
    return result;
  }
  
  /// Match the arithmetic precedence levels shared by the script grammar
  int arithmeticPrecedence(const std::string_view token)
  {
    int precedence=0;
    if(token=="+" or token=="-")
      precedence=1;
    else if(token=="*" or token=="/" or token=="%")
      precedence=2;
    
    return precedence;
  }
  
  /// Select the same registered overload as runtime dispatch, without invoking it
  template <typename Operator>
  CompletionBinding binaryResult(const CompletionBinding& left,
                                const CompletionBinding& right)
  {
    CompletionBinding result;
    if(left.cppType and right.cppType)
      {
        const auto host=left.hostType?left.hostType:right.hostType;
        if(host)
          {
            if(const auto operation=host->binaryOperations.find(typeid(Operator));operation!=host->binaryOperations.end())
              {
                const auto& overloads=left.hostType?operation->second.rightOperands:operation->second.leftOperands;
                const auto found=overloads.find(left.hostType?*right.cppType:*left.cppType);
                if(found!=overloads.end())
                  result=resultBinding(found->second.resultType);
              }
          }
        else
          {
            const auto numeric=[](const std::type_index type)
            {
              return type==typeid(int) or type==typeid(int64_t) or type==typeid(double);
            };
            if(numeric(*left.cppType) and numeric(*right.cppType))
              {
                const bool real=*left.cppType==typeid(double) or *right.cppType==typeid(double);
                if(not real or not std::is_same_v<Operator,RemainderOperator>)
                  result=resultBinding({real?typeid(double):
                    (*left.cppType==typeid(int64_t) or *right.cppType==typeid(int64_t)?typeid(int64_t):typeid(int))});
              }
          }
      }
    
    return result;
  }
  
  /// Small tolerant type reader for calls, grouping, and registered arithmetic expressions
  struct CompletionExpression
  {
    /// Visible binding metadata borrowed for the duration of this query
    const std::map<std::string,CompletionBinding>& bindings;
    
    /// Complete initializer tokens; no parser or semantic actions are invoked
    std::span<const std::string_view> tokens;
    
    /// Next unconsumed token
    size_t position{0};
    
    /// Reject malformed expressions and cap recursion in untrusted editor input
    bool valid{true};
    
    /// Read a primary, its calls and member types, or a numeric unary operation
    CompletionBinding primary(const size_t nesting)
    {
      CompletionBinding value;
      if(position==tokens.size() or nesting>128)
        valid=false;
      else
        {
          const auto token=tokens[position++];
          if(token=="(")
            {
              value=expression(1,nesting+1);
              if(position<tokens.size() and tokens[position]==")")
                position++;
              else
                valid=false;
            }
          else if(token=="+" or token=="-")
            {
              value=primary(nesting+1);
              if(value.hostType or not value.cppType or
                 (*value.cppType!=typeid(int) and *value.cppType!=typeid(int64_t) and *value.cppType!=typeid(double)))
                value={};
            }
          else if(token.starts_with('"'))
            value=resultBinding({typeid(std::string)});
          else if(token.front()>='0' and token.front()<='9')
            {
              int64_t integer{};
              const auto parsed=std::from_chars(token.data(),token.data()+token.size(),integer);
              if(parsed.ec==std::errc{} and parsed.ptr==token.data()+token.size())
                value=resultBinding({integer<=std::numeric_limits<int>::max()?typeid(int):typeid(int64_t)});
              else if(token.find_first_of(".eE")!=std::string_view::npos)
                value=resultBinding({typeid(double)});
            }
          else if(const auto found=bindings.find(std::string(token));found!=bindings.end())
            value=found->second;
          
          while(valid and position<tokens.size() and (tokens[position]=="(" or tokens[position]=="."))
            {
              if(tokens[position]==".")
                {
                  position++;
                  CompletionBinding member;
                  if(position<tokens.size() and identifier(tokens[position]))
                    {
                      const std::string name(tokens[position++]);
                      if(value.hostType)
                        {
                          const auto found=value.hostType->memberAccessorsProvider.find(name);
                          if(found!=value.hostType->memberAccessorsProvider.end())
                            member=resultBinding(found->second.valueType);
                        }
                      else if(const auto found=value.members.find(name);found!=value.members.end())
                        member.info=found->second;
                    }
                  else
                    valid=false;
                  
                  value=std::move(member);
                }
              else
                {
                  // Fixed native return types do not require evaluating or typing individual arguments
                  std::vector<std::string_view> closing;
                  do
                    {
                      const auto argument=tokens[position++];
                      if(argument=="(" or argument=="[")
                        closing.push_back(argument=="("?")":"]");
                      else if(argument==")" or argument=="]")
                        {
                          valid=not closing.empty() and closing.back()==argument;
                          if(valid)
                            closing.pop_back();
                        }
                    }
                  while(valid and not closing.empty() and position<tokens.size());
                  
                  valid=valid and closing.empty();
                  if(value.scriptFunction)
                    value=scriptReturn(*value.scriptFunction,bindings);
                  else
                    value=value.callResult?*value.callResult:CompletionBinding{};
                }
            }
        }
      
      return value;
    }
    
    /// Infer left-associative arithmetic with the grammar's precedence, without computing values
    CompletionBinding expression(const int minimumPrecedence,
                                 const size_t nesting)
    {
      auto left=primary(nesting);
      while(valid and position<tokens.size() and arithmeticPrecedence(tokens[position])>=minimumPrecedence)
        {
          const auto operation=tokens[position++];
          const auto right=expression(arithmeticPrecedence(operation)+1,nesting+1);
          if(operation=="+")
            left=binaryResult<SumOperator>(left,right);
          else if(operation=="-")
            left=binaryResult<DifferenceOperator>(left,right);
          else if(operation=="*")
            left=binaryResult<ProductOperator>(left,right);
          else if(operation=="/")
            left=binaryResult<DivisionOperator>(left,right);
          else
            left=binaryResult<RemainderOperator>(left,right);
        }
      
      return left;
    }
  };
  
  /// Infer a complete initializer; unknown or malformed expressions retain no guessed result
  CompletionBinding inferBinding(const std::map<std::string,CompletionBinding>& bindings,
                                 const std::span<const std::string_view> expression)
  {
    CompletionExpression reader{bindings,expression};
    auto result=reader.expression(1,0);
    if(not reader.valid or reader.position!=expression.size())
      result={};
    
    return result;
  }
  
  /// Compare known value types and preserve callable identity when merging local bindings
  bool sameCompletionType(const CompletionBinding& left,
                          const CompletionBinding& right)
  {
    return left.cppType and left.cppType==right.cppType and left.hostType==right.hostType and
      left.callResult==right.callResult and left.scriptFunction==right.scriptFunction;
  }
  
  /// Tolerant statement reader tracking return types and the remaining fallthrough path
  struct CompletionBody
  {
    /// Body tokens borrowed for this completion query
    std::span<const std::string_view> tokens;
    
    /// Next token to consume
    size_t position{0};
    
    /// Unsupported reachable syntax makes the whole result unknown
    bool valid{true};
    
    /// Common type of all reachable returns encountered so far
    std::optional<CompletionBinding> returned;
    
    /// Add a return without treating a missing or conflicting type as agreement
    void addReturn(const CompletionBinding& value)
    {
      if(not returned)
        returned=value;
      else if(not sameCompletionType(*returned,value))
        returned=CompletionBinding{};
    }
    
    /// Read statements until the current block ends, ignoring unreachable return values
    bool sequence(std::map<std::string,CompletionBinding>& locals,
                  const bool reachable,
                  const size_t nesting)
    {
      bool fallsThrough=reachable;
      while(valid and position<tokens.size() and tokens[position]!="}")
        fallsThrough=statement(locals,fallsThrough,nesting+1);
      
      return fallsThrough;
    }
    
    /// Read a block, conditional, return, or simple assignment without evaluating expressions
    bool statement(std::map<std::string,CompletionBinding>& locals,
                   const bool reachable,
                   const size_t nesting)
    {
      bool fallsThrough=reachable;
      if(position==tokens.size() or nesting>128)
        valid=false;
      else if(tokens[position]=="{")
        {
          position++;
          auto inner=locals;
          fallsThrough=sequence(inner,reachable,nesting);
          if(position<tokens.size() and tokens[position]=="}")
            position++;
          else
            valid=false;
          
          for(auto& [name,value] : locals)
            value=inner.at(name);
        }
      else if(tokens[position]=="if")
        {
          position++;
          if(position==tokens.size() or tokens[position]!="(")
            valid=false;
          else
            {
              size_t depth=0;
              do
                {
                  const auto token=tokens[position++];
                  if(token=="(")
                    depth++;
                  else if(token==")")
                    depth--;
                  else if(token=="=" or token=="+=" or token=="-=" or token=="*=" or token=="/=" or
                          token=="++" or token=="--" or
                          (identifier(token) and position<tokens.size() and tokens[position]=="("))
                    valid=false;
                }
              while(valid and depth>0 and position<tokens.size());
              valid=valid and depth==0;
            }
          
          if(valid)
            {
              auto left=locals;
              auto right=locals;
              const bool leftFalls=statement(left,reachable,nesting+1);
              bool rightFalls=reachable;
              if(valid and position<tokens.size() and tokens[position]=="else")
                {
                  position++;
                  rightFalls=statement(right,reachable,nesting+1);
                }
              
              // Only paths reaching the next statement contribute its local types
              for(auto& [name,value] : locals)
                if(leftFalls and rightFalls)
                  value=sameCompletionType(left.at(name),right.at(name))?left.at(name):CompletionBinding{};
                else if(leftFalls)
                  value=left.at(name);
                else if(rightFalls)
                  value=right.at(name);
              
              fallsThrough=leftFalls or rightFalls;
            }
        }
      else
        {
          const size_t begin=position;
          size_t depth=0;
          while(valid and position<tokens.size() and not (tokens[position]==";" and depth==0))
            {
              const auto token=tokens[position++];
              if(token=="(" or token=="[")
                depth++;
              else if(token==")" or token=="]")
                {
                  valid=depth>0;
                  if(valid)
                    depth--;
                }
              else if(token=="{" or token=="}")
                valid=false;
            }
          
          if(valid and position<tokens.size())
            {
              const auto expression=tokens.subspan(begin,position-begin);
              position++;
              if(reachable and not expression.empty())
                {
                  if(expression.front()=="return")
                    {
                      addReturn(inferBinding(locals,expression.subspan(1)));
                      fallsThrough=false;
                    }
                  else if(expression.size()>2 and identifier(expression[0]) and expression[1]=="=")
                    locals[std::string(expression[0])]=inferBinding(locals,expression.subspan(2));
                  else
                    valid=false;
                }
            }
          else
            valid=false;
        }
      
      return fallsThrough;
    }
  };
  
  CompletionBinding scriptReturn(const CompletionFunction& function,
                                 const std::map<std::string,CompletionBinding>& bindings)
  {
    std::map<std::string,CompletionBinding> locals;
    for(const auto& [name,result] : function.hostResults)
      if(const auto found=bindings.find(name);found!=bindings.end() and found->second.callResult==result)
        locals.emplace(name,found->second);
    
    CompletionBody reader{function.body};
    const bool fallsThrough=reader.sequence(locals,true,0);
    CompletionBinding result;
    if(reader.valid and reader.position==function.body.size() and not fallsThrough and reader.returned)
      result=*reader.returned;
    
    return result;
  }
  
  /// Record completed assignments, forgetting type information for unsupported expressions
  void recordAssignments(CompletionScope& scope,
                         const std::span<const std::string_view> statement,
                         const ScriptSource& source)
  {
    const bool conditional=not statement.empty() and
      (statement.front()=="if" or statement.front()=="for" or statement.front()=="else");
    size_t depth=0;
    for(size_t i=0;i<statement.size();i++)
      {
        if(statement[i]=="(" or statement[i]=="[")
          depth++;
        else if((statement[i]==")" or statement[i]=="]") and depth>0)
          depth--;
        else if(i>0 and identifier(statement[i-1]) and
                (i<2 or statement[i-2]!=".") and
                (statement[i]=="=" or statement[i]=="+=" or statement[i]=="-=" or
                 statement[i]=="*=" or statement[i]=="/=") and
                (not conditional or scope.bindings.contains(std::string(statement[i-1]))))
          {
            const std::string name(statement[i-1]);
            CompletionBinding binding;
            if(depth==0 and i==1 and statement[i]=="=")
              binding=inferBinding(scope.bindings,statement.subspan(i+1));
            
            const auto location=source.locate(SourceCursor{
              .offset=static_cast<size_t>(statement[i-1].data()-source.text.data())});
            if(not binding.info.callable or not binding.info.location)
              binding.info.location=ScriptLocation{std::string(location.source),location.cursor};
            
            scope.bindings.insert_or_assign(name,std::move(binding));
            scope.written.insert(name);
          }
      }
  }
  
  /// Reconstruct straight-line bindings and brace visibility without invoking the parser or evaluator
  std::map<std::string,CompletionBinding> completionBindings(const Environment& environment,
                                                            const std::vector<std::string_view>& tokens,
                                                            const ScriptSource& source)
  {
    std::vector<CompletionScope> scopes(1);
    for(const auto& name : environment.entityNames())
      {
        CompletionBinding binding{.info=environment.entityInfo(name),.documentation=formatHelp(environment,name)};
        const Environment* scope=&environment;
        while(scope and not scope->varTable.contains(name))
          scope=scope->parent.get();
        
        if(const auto& value=scope->varTable.at(name);value)
          {
            binding.cppType=std::visit([](const auto& payload)
            {
              return std::type_index(typeid(payload));
            },value->data);
            
            if(const auto* object=std::get_if<HObject>(&value->data);object and object->info)
              {
                binding.hostType=object->info;
                binding.cppType=object->info->cppType;
                binding.members=hostMembers(*object->info);
              }
            else if(const auto* object=std::get_if<Struct>(&value->data);object and object->structDefNode)
              binding.members=scriptMembers(*object->structDefNode);
            else if(const auto* function=std::get_if<HostFunction>(&value->data);function and function->resultType.cppType)
              binding.callResult=std::make_shared<CompletionBinding>(resultBinding(function->resultType));
            else if(const auto* definition=std::get_if<StructDef>(&value->data);definition and definition->structDefNode)
              {
                binding.callResult=std::make_shared<CompletionBinding>();
                binding.callResult->members=scriptMembers(*definition->structDefNode);
                binding.callResult->info.kind="script struct instance '"+definition->structDefNode->name+"'";
              }
          }
        
        scopes.back().bindings.emplace(name,std::move(binding));
      }
    
    size_t depth=0;
    const std::span<const std::string_view> all(tokens);
    for(size_t i=0;i<tokens.size();i++)
      {
        if(tokens[i]=="{")
          {
            const auto header=all.subspan(scopes.back().begin,i-scopes.back().begin);
            const bool deferred=std::find(header.begin(),header.end(),"fun")!=header.end() or
              std::find(header.begin(),header.end(),"lambda")!=header.end() or
              std::find(header.begin(),header.end(),"struct")!=header.end();
            // Record the assigned name before scanning a deferred lambda body
            if(deferred and header.size()>1 and header[1]=="=")
              recordAssignments(scopes.back(),header,source);
            
            CompletionScope child{.bindings=scopes.back().bindings,.begin=i+1,.parentDepth=depth,
              .kind=header.empty()?CompletionScope::Kind::Plain:
              (deferred?CompletionScope::Kind::Deferred:CompletionScope::Kind::Conditional)};
            
            if(not deferred and not header.empty())
              recordAssignments(child,header.subspan(1),source);
            
            if(deferred)
              {
                for(size_t j=0;j<header.size();j++)
                  if(((header[j]=="fun" or header[j]=="struct") and j==0 and depth==0) or header[j]=="lambda")
                    {
                      const bool structure=header[j]=="struct";
                      const bool lambda=header[j]=="lambda";
                      size_t open=j+1;
                      std::string name;
                      if(not lambda and j==0 and open<header.size() and identifier(header[open]))
                        name=header[open++];
                      else if(lambda and j==2 and header[1]=="=" and identifier(header[0]))
                        name=header[0];
                      
                      CompletionBinding definition;
                      definition.info={.kind=structure?"script struct constructor":
                        (lambda?"script lambda":"embedded function '"+name+"'"),
                        .description=completionDescription(source,all,scopes.back().begin+j),
                        .callable=true,.argumentForms=FunctionArgumentForms(1),
                        .location=completionLocation(source,header[j])};
                      if(lambda and j==2 and definition.info.description.empty())
                        definition.info.description=completionDescription(source,all,scopes.back().begin);
                      
                      std::set<std::string> parameterNames;
                      child.definitionName=name;
                      child.structDefinition=structure;
                      child.lambdaDefinition=lambda;
                      if(structure)
                        {
                          definition.callResult=std::make_shared<CompletionBinding>();
                          definition.callResult->info.kind="script struct instance '"+name+"'";
                        }
                      else if(open<header.size() and header[open]=="(" and header.back()==")")
                        for(const auto& parameter : completionParameters(header.subspan(open+1,header.size()-open-2),source))
                          {
                            definition.info.argumentForms[0].push_back(parameter.form);
                            parameterNames.insert(parameter.name);
                            child.bindings[parameter.name]=CompletionBinding{
                              .info={.kind="function parameter",.location=parameter.location}};
                          }
                      
                      if(not structure)
                        {
                          definition.scriptFunction=std::make_shared<CompletionFunction>();
                          definition.scriptFunction->body=all.subspan(i+1);
                          for(const auto& [visible,binding] : child.bindings)
                            if(binding.callResult and not binding.info.location)
                              definition.scriptFunction->hostResults.emplace(visible,binding.callResult);
                        }
                      
                      definition.documentation=completionDocumentation(name.empty()?"lambda":name,definition.info);
                      if(not name.empty())
                        {
                          scopes.back().bindings[name]=definition;
                          scopes.back().written.insert(name);
                          if(not parameterNames.contains(name))
                            child.bindings[name]=definition;
                        }
                      child.definition=std::move(definition);
                    }
              }
            
            scopes.push_back(std::move(child));
            depth=0;
          }
        else if(tokens[i]=="}" and scopes.size()>1)
          {
            auto child=std::move(scopes.back());
            scopes.pop_back();
            auto& parent=scopes.back();
            if(child.definition and not child.definitionName.empty())
              {
                if(child.definition->scriptFunction)
                  {
                    auto& body=child.definition->scriptFunction->body;
                    body=body.first(all.data()+i-body.data());
                  }
                
                if(child.lambdaDefinition and i+1<tokens.size() and tokens[i+1]!=";")
                  child.definition=CompletionBinding{};
                
                child.definition->documentation=completionDocumentation(child.definitionName,child.definition->info);
                parent.bindings[child.definitionName]=std::move(*child.definition);
                parent.written.insert(child.definitionName);
              }
            
            if(child.kind!=CompletionScope::Kind::Deferred)
              for(const auto& name : child.written)
                if(parent.bindings.contains(name))
                  {
                    parent.bindings[name]=child.kind==CompletionScope::Kind::Plain?
                      child.bindings.at(name):CompletionBinding{};
                    parent.written.insert(name);
                  }
            
            parent.begin=i+1;
            depth=child.parentDepth;
          }
        else if(tokens[i]=="(" or tokens[i]=="[")
          depth++;
        else if((tokens[i]==")" or tokens[i]=="]") and depth>0)
          depth--;
        else if(tokens[i]==";" and depth==0)
          {
            auto& scope=scopes.back();
            const auto statement=all.subspan(scope.begin,i-scope.begin);
            if(scope.structDefinition and scope.definition)
              {
                if(not statement.empty() and identifier(statement.front()) and
                   (statement.size()==1 or (statement.size()>2 and statement[1]=="=")))
                  {
                    const std::string name(statement.front());
                    EntityInfo member{.kind="script member",
                      .description=completionDescription(source,all,scope.begin),
                      .location=completionLocation(source,statement.front())};
                    scope.definition->callResult->members.insert_or_assign(name,std::move(member));
                    std::string form=name;
                    if(statement.size()>2)
                      {
                        const auto first=statement[2];
                        const auto last=statement.back();
                        form="["+name+" = "+std::string(first.data(),last.data()+last.size()-first.data())+"]";
                      }
                    scope.definition->info.argumentForms[0].push_back(std::move(form));
                  }
              }
            else
              recordAssignments(scope,statement,source);
            
            scope.begin=i+1;
          }
      }
    
    return std::move(scopes.back().bindings);
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
      ScriptSource source;
      source.append(text,"<script>");
      result=completeScriptAt(environment,source,byteOffset);
    }
  
  return result;
}

ScriptCompletion completeScriptAt(const Environment& environment,
                                  const ScriptSource& source,
                                  const size_t byteOffset)
{
  ScriptCompletion result;
  const std::string_view text=source.text;
  if(byteOffset<=text.size())
    {
      const std::string_view scanned=text;
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
                {
                  state=State::Code;
                  const size_t begin=tokens.back().data()-scanned.data();
                  tokens.back()=scanned.substr(begin,i+1-begin);
                }
              
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
              tokens.push_back(scanned.substr(i,1));
            }
          else if(c>='0' and c<='9')
            {
              double number{};
              const auto parsed=std::from_chars(text.data()+i,text.data()+result.begin,number);
              const size_t length=parsed.ptr-(text.data()+i);
              tokens.push_back(scanned.substr(i,length));
              i+=length-1;
            }
          else if(identifierByte(c))
            {
              const size_t begin=i;
              while(i+1<result.begin and identifierByte(text[i+1]))
                i++;
              
              tokens.push_back(scanned.substr(begin,i+1-begin));
            }
          else if(c!=' ' and c!='\t' and c!='\r' and c!='\n' and c!='\f' and c!='\v')
            {
              if(c=='(' or c=='[' or c=='{')
                delimiters.push_back(tokens.size());
              else if((c==')' or c==']' or c=='}') and not delimiters.empty())
                delimiters.pop_back();
              
              const size_t length=i+1<result.begin and
                ((next=='=' and (c=='=' or c=='!' or c=='<' or c=='>' or c=='+' or c=='-' or c=='*' or c=='/')) or
                 (next==c and (c=='+' or c=='-' or c=='&' or c=='|')))?2:1;
              tokens.push_back(scanned.substr(i,length));
              i+=length-1;
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
          const auto bindings=completionBindings(environment,tokens,source);
          const auto infoOf=[&](const std::string_view name)->std::optional<EntityInfo>
          {
            std::optional<EntityInfo> info;
            if(const auto found=bindings.find(std::string(name));found!=bindings.end())
              info=found->second.info;
            
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
                                  add(argument,EntityInfo{"parameter"},bindings.at(std::string(tokens[open-1])).documentation);
                              }
                    }
                }
              else if(tokens.size()>1)
                {
                  // Find the receiver expression within its statement or enclosing call argument
                  const size_t end=tokens.size()-1;
                  size_t begin=end;
                  size_t depth=0;
                  bool searching=true;
                  while(begin>0 and searching)
                    {
                      const auto token=tokens[begin-1];
                      if(token==")" or token=="]")
                        depth++;
                      else if(token=="(" or token=="[")
                        {
                          if(depth>0)
                            depth--;
                          else
                            searching=false;
                        }
                      else if(depth==0 and (token==";" or token=="{" or token=="}" or
                                              token=="=" or token=="," or token=="return"))
                        searching=false;
                      
                      if(searching)
                        begin--;
                    }
                  
                  const auto object=inferBinding(bindings,std::span<const std::string_view>(tokens).subspan(begin,end-begin));
                  for(const auto& [name,info] : object.members)
                    add(name,info,info.description);
                }
            }
          else
            {
              for(const auto& [name,binding] : bindings)
                add(name,binding.info,binding.documentation);
              
              for(const std::string keyword : {"fun","lambda","struct","if","else","for","return","and","or"})
                if(not bindings.contains(keyword))
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
