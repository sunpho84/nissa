#include "scriptCompletion.hpp"
#include "scriptAnalysis.hpp"

#include <sstream>
#include <map>

namespace
{
  /// Static value metadata; no runtime payload or getter is accessed
  struct CompletionBinding
  {
    /// Classification and declaration position
    EntityInfo info{"script variable (unknown type)"};
    
    /// Complete help supplied to editor clients
    std::string documentation;
    
    /// Statically available members
    std::map<std::string,EntityInfo> members;
    
    /// Native scalar or payload identity
    std::optional<std::type_index> cppType{};
    
    /// Shared host type metadata
    std::shared_ptr<const HObjectInfo> hostType;
    
    /// Fixed result metadata for native functions and constructors
    std::shared_ptr<CompletionBinding> callResult;
    
    /// Original AST body for a script callable
    std::shared_ptr<const FuncNode> scriptFunction;
    
    /// Script constructor or instance identity, borrowed from this revision's AST
    const StructDefNode* scriptStruct{nullptr};
    
    /// Inferred member types for this particular constructed value
    std::map<std::string,std::shared_ptr<CompletionBinding>> memberTypes;
    
    /// Callable bindings visible at definition, checked again when inferring a call
    std::map<std::string,std::shared_ptr<CompletionBinding>> definitionBindings;
  };
  
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
  
  /// Compare static value types, including the identity of callable values
  bool sameBinding(const CompletionBinding& left,
                   const CompletionBinding& right,
                   size_t* budget=nullptr)
  {
    size_t remaining=4096;
    if(not budget)
      budget=&remaining;
    bool same=*budget>0 and ((left.cppType and left.cppType==right.cppType) or
               (left.callResult and left.callResult==right.callResult) or
               (left.scriptFunction and left.scriptFunction==right.scriptFunction) or
               (left.scriptStruct and left.scriptStruct==right.scriptStruct)) and
      left.hostType==right.hostType and left.callResult==right.callResult and
      left.scriptFunction==right.scriptFunction and left.scriptStruct==right.scriptStruct and
      left.memberTypes.size()==right.memberTypes.size();
    if(*budget>0)
      (*budget)--;
    for(const auto& [name,type] : left.memberTypes)
      if(const auto other=right.memberTypes.find(name);same and
         (other==right.memberTypes.end() or
          (type!=other->second and not sameBinding(*type,*other->second,budget))))
        same=false;
    
    return same;
  }
  
  /// Static AST visitor shared by scope reconstruction and script-return inference
  struct CompletionAnalyzer
  {
    /// Owned AST and source ranges for this editor revision
    const ParsedScript& script;
    
    /// Current lexical scope, independently of runtime storage
    std::map<std::string,CompletionBinding> bindings;
    
    /// Cursor target in the retained source, independent of AST mutation
    std::optional<SourceRange> target{};
    
    /// Range of the node currently being visited
    std::optional<SourceRange> activeRange{};
    
    /// Check adjacency using source whitespace, without recognizing syntax again
    bool whitespace(const size_t begin,
                    const size_t end) const
    {
      return begin<=end and end<=script.source.text.size() and
        std::string_view(script.source.text).substr(begin,end-begin).find_first_not_of(" \t\r\n\f\v")==std::string_view::npos;
    }
    
    /// Whether the current node contains a gap cursor
    bool gapHere() const
    {
      return not returnMode and target and target->begin==target->end and activeRange and
        activeRange->begin.offset<=target->begin.offset and target->end.offset<=activeRange->end.offset;
    }
    
    /// Match actual identifier positions, or a recovery placeholder adjacent to a gap
    bool matches(const std::string& name,
                 const std::optional<ScriptLocation>& location={}) const
    {
      bool result=false;
      if(not returnMode and target)
        {
          if(script.completion)
            result=name==completionMarker;
          else if(location and activeRange and
                  (gapHere() or (activeRange->begin.offset<=target->begin.offset and
                                 target->end.offset<=activeRange->end.offset) or
                   (name==recoveryMarker and whitespace(target->end.offset,activeRange->begin.offset))))
            {
              const auto position=script.source.locate(target->begin);
              result=location->source==position.source and location->cursor.offset==position.cursor.offset;
              if(name==recoveryMarker and target->begin==target->end)
                {
                  // Synthetic names may follow the final newline added when pieces are joined
                  const auto end=script.source.locate(activeRange->end);
                  result=location->source==end.source and location->cursor.offset==end.cursor.offset and
                    whitespace(target->begin.offset,activeRange->end.offset);
                }
            }
        }
      
      return result;
    }
    
    /// Marker was reached; bindings now describe its scope
    bool found{false};
    
    /// Member or named-argument completion, rather than lexical names
    bool restricted{false};
    
    /// Exact members or arguments offered at the marker
    std::map<std::string,CompletionBinding> choices;
    
    /// Function-return analysis never captures a marker from a deferred body
    bool returnMode{false};
    
    /// A reachable path still continues to the next statement
    bool fallsThrough{true};
    
    /// Unsupported effects make return inference unknown
    bool invalid{false};
    
    /// Common known type of all reachable return expressions
    std::optional<CompletionBinding> returned;
    
    /// Bound recursive AST walks and return inference
    size_t depth{0};
    
    /// Enclosing names whose assignment would mutate a captured binding
    std::set<std::string> capturedNames;
    
    /// Captured bindings assigned by this body or a nested helper
    std::set<std::string> capturedWrites;
    
    /// Shared work limit across call specializations and branches, protecting editor responsiveness
    std::shared_ptr<size_t> remaining{std::make_shared<size_t>(16384)};
    
    /// Read a node using static metadata only
    CompletionBinding visit(const std::shared_ptr<ASTNode>& node)
    {
      CompletionBinding result;
      if(node and *remaining==0)
        invalid=true;
      if(node and not found and depth<128 and *remaining>0 and (not returnMode or fallsThrough))
        {
          const auto previous=activeRange;
          const auto range=script.ranges.find(node.get());
          activeRange=range==script.ranges.end()?std::optional<SourceRange>{}:range->second;
          if(not returnMode and not script.completion and target and target->begin==target->end and activeRange and
             whitespace(target->begin.offset,activeRange->begin.offset))
            found=true;
          depth++;
          (*remaining)--;
          if(not found)
            result=std::visit([this](const auto& value)
          {
            return analyze(value);
          },*node);
          depth--;
          const bool closedBlock=std::holds_alternative<ASTNodesNode>(*node) and activeRange and target and
            target->begin.offset==activeRange->end.offset and target->begin.offset>0 and
            script.source.text[target->begin.offset-1]=='}';
          if(not script.completion and gapHere() and not closedBlock)
            found=true;
          activeRange=previous;
        }
      
      return result;
    }
    
    /// Require an overload when the AST gains a new node alternative
    template <typename T>
    CompletionBinding analyze(const T&)=delete;
    
    /// Resolve a name without creating or reading its runtime value
    CompletionBinding analyze(const IdNode& node)
    {
      CompletionBinding result;
      if(matches(node.name,node.location))
        found=true;
      else if(const auto entry=bindings.find(node.name);entry!=bindings.end())
        result=entry->second;
      
      return result;
    }
    
    /// Literal types require no value evaluation
    CompletionBinding analyze(const ValueNode& node)
    {
      return std::visit([](const auto& value)
      {
        return resultBinding({typeid(value)});
      },node.value.data);
    }
    
    /// Traverse a statement list, optionally preserving lexical block boundaries
    void statements(const ASTNodesNode& node,
                    const bool block)
    {
      auto outer=bindings;
      for(const auto& statement : node.list)
        {
          if(returnMode and fallsThrough and not std::holds_alternative<AssignNode>(*statement) and
             not std::holds_alternative<ReturnNode>(*statement) and not std::holds_alternative<IfNode>(*statement) and
             not std::holds_alternative<ASTNodesNode>(*statement) and not std::holds_alternative<ValueNode>(*statement))
            invalid=true;
          visit(statement);
        }
      if(not script.completion and gapHere() and
         (not block or target->begin.offset<activeRange->end.offset))
        found=true;
      if(block and not found)
        {
          for(auto& [name,value] : outer)
            value=bindings.at(name);
          bindings=std::move(outer);
        }
    }
    
    /// Braces introduce a child scope, just as in the runtime evaluator
    CompletionBinding analyze(const ASTNodesNode& node)
    {
      statements(node,true);
      
      return {};
    }
    
    /// Retain assignment types and callable aliases, without storing script values
    CompletionBinding analyze(const AssignNode& node)
    {
      auto value=visit(node.rhs);
      const auto* name=std::get_if<IdNode>(node.lhs.get());
      if(name and matches(name->name,name->location))
        found=true;
      else if(not found and name and name->name!=recoveryMarker)
        {
          if(returnMode and capturedNames.contains(name->name))
            {
              invalid=true;
              capturedWrites.insert(name->name);
            }
          if(node.mode!=AssignMode::PLAIN)
            value={};
          if(not value.info.callable or not value.info.location)
            value.info.location=name->location;
          bindings[name->name]=value;
        }
      else if(not name)
        {
          visit(node.lhs);
          if(std::holds_alternative<MemberAccessNode>(*node.lhs))
            {
              // Struct copies can share member cells; discard all potentially aliased shapes
              for(auto& [visible,binding] : bindings)
                binding.memberTypes.clear();
              if(returnMode)
                invalid=true;
            }
        }
      
      return value;
    }
    
    /// Registered binary operators expose the same type-dispatch metadata as execution
    CompletionBinding analyze(const BinOpNode& node)
    {
      const auto left=visit(node.arg1);
      const auto right=visit(node.arg2);
      CompletionBinding result;
      if(node.operatorType==binaryOperatorType<SumOperator>())
        result=binaryResult<SumOperator>(left,right);
      else if(node.operatorType==binaryOperatorType<DifferenceOperator>())
        result=binaryResult<DifferenceOperator>(left,right);
      else if(node.operatorType==binaryOperatorType<ProductOperator>())
        result=binaryResult<ProductOperator>(left,right);
      else if(node.operatorType==binaryOperatorType<DivisionOperator>())
        result=binaryResult<DivisionOperator>(left,right);
      else if(node.operatorType==binaryOperatorType<RemainderOperator>())
        result=binaryResult<RemainderOperator>(left,right);
      
      return result;
    }
    
    /// Numeric unary operations preserve the argument type; logical negation returns an integer
    CompletionBinding analyze(const UnOpNode& node)
    {
      auto value=visit(node.arg);
      if(node.kind==UnOpNode::Kind::Not)
        value=resultBinding({typeid(int)});
      else if(value.hostType)
        value={};
      
      return value;
    }
    
    /// A reference value is not itself its referent; only explicit reference arguments specialize
    CompletionBinding analyze(const UnaryRefNode& node)
    {
      visit(node.arg);
      
      return {};
    }
    
    /// Mutation invalidates a known local type
    CompletionBinding analyze(const PrePostfixOpNode& node)
    {
      visit(node.arg);
      if(const auto* id=std::get_if<IdNode>(node.arg.get());id and bindings.contains(id->name))
        bindings[id->name]={};
      else if(std::holds_alternative<MemberAccessNode>(*node.arg))
        {
          for(auto& [name,binding] : bindings)
            binding.memberTypes.clear();
          if(returnMode)
            invalid=true;
        }
      
      return {};
    }
    
    /// Member lookup never invokes a registered accessor
    CompletionBinding analyze(const MemberAccessNode& node)
    {
      const auto base=visit(node.base);
      CompletionBinding result;
      if(matches(node.member,node.memberLocation))
        {
          found=true;
          restricted=true;
          for(const auto& [name,info] : base.members)
            choices.emplace(name,CompletionBinding{.info=info,.documentation=info.description});
        }
      else if(base.hostType)
        {
          if(const auto member=base.hostType->memberAccessorsProvider.find(node.member);
             member!=base.hostType->memberAccessorsProvider.end())
            result=resultBinding(member->second.valueType);
        }
      else if(const auto member=base.memberTypes.find(node.member);member!=base.memberTypes.end())
        result=*member->second;
      else if(const auto member=base.members.find(node.member);member!=base.members.end())
        result.info=member->second;
      
      return result;
    }
    
    /// Indexed container result types remain unknown
    CompletionBinding analyze(const SubscribeNode& node)
    {
      visit(node.base);
      visit(node.subscr);
      
      return {};
    }
    
    /// Recover source spelling for unevaluated default expressions
    std::string spelling(const std::shared_ptr<ASTNode>& node) const
    {
      std::string result="<expression>";
      if(const auto range=script.ranges.find(node.get());range!=script.ranges.end())
        result=script.source.text.substr(range->second.begin.offset,
                                         range->second.end.offset-range->second.begin.offset);
      
      return result;
    }
    
    /// Build a script callable from its original AST and retain visible callable metadata
    CompletionBinding function(const std::shared_ptr<const FuncNode>& node,
                               const std::string& name)
    {
      CompletionBinding result;
      result.info={.kind=name.empty()?"script lambda":"embedded function '"+name+"'",
        .description=node->description,.callable=true,.argumentForms=FunctionArgumentForms(1),.location=node->location};
      for(const auto& parameter : node->pars.list)
        {
          std::string form=parameter.isRef?"&":"";
          form+=parameter.name==FuncParNode::VARIADIC_PAR_NAME?"...":parameter.name;
          if(parameter.def and parameter.name!=FuncParNode::VARIADIC_PAR_NAME)
            form="["+form+" = "+spelling(parameter.def)+"]";
          result.info.argumentForms[0].push_back(std::move(form));
        }
      result.documentation=completionDocumentation(name.empty()?"lambda":name,result.info);
      result.scriptFunction=node;
      for(const auto& [visible,binding] : bindings)
        if(binding.callResult or binding.scriptFunction)
          result.definitionBindings.emplace(visible,std::make_shared<CompletionBinding>(binding));
      for(const auto& parameter : node->pars.list)
        result.definitionBindings.erase(parameter.name);
      
      return result;
    }
    
    /// Inspect a deferred body for the editor marker without affecting the enclosing scope
    void inspectFunction(const FuncNode& node)
    {
      if(not returnMode)
        {
          auto body=*this;
          for(const auto& parameter : node.pars.list)
            {
              body.visit(parameter.def);
              if(not body.found)
                body.bindings[parameter.name]=CompletionBinding{
                  .info={.kind="function parameter",.location=parameter.location}};
            }
          body.visit(node.body);
          if(body.found)
            adopt(body);
        }
    }
    
    /// Copy the scope and completion target reached by a nested visitor
    void adopt(const CompletionAnalyzer& child)
    {
      bindings=child.bindings;
      found=child.found;
      restricted=child.restricted;
      choices=child.choices;
    }
    
    /// Named declarations install the callable before inspecting its body
    CompletionBinding analyze(const FuncDefNode& node)
    {
      auto value=function(node.fun,node.name);
      if(matches(node.name,node.location))
        found=true;
      else if(node.name!=recoveryMarker)
        bindings[node.name]=value;
      inspectFunction(*node.fun);
      
      return value;
    }
    
    /// Anonymous functions share the same declaration and return analysis
    CompletionBinding analyze(const FuncNode& node)
    {
      auto value=function(std::make_shared<FuncNode>(node),"");
      inspectFunction(node);
      
      return value;
    }
    
    /// Argument metadata, computed once in source order through the normal AST visitor
    struct Argument
    {
      /// Empty for a positional argument
      std::string name;
      
      /// Metadata of the expression, without evaluating it
      CompletionBinding type;
      
      /// Whether the AST explicitly passes a script reference
      bool reference{false};
    };
    
    /// Match script arguments in runtime order: positional first, then named, then defaults
    template <typename Parameter>
    std::optional<std::vector<CompletionBinding>> arguments(const std::vector<Parameter>& parameters,
                                                          const bool variadic,
                                                          const std::vector<Argument>& supplied)
    {
      std::vector<std::optional<CompletionBinding>> slots(parameters.size());
      bool valid=true;
      if(variadic)
        slots.back()=resultBinding({typeid(ValuesList)});
      for(const bool named : {false,true})
        for(const auto& argument : supplied)
          if(named==not argument.name.empty())
            {
              size_t index=0;
              if(named)
                while(index<parameters.size() and parameters[index].name!=argument.name)
                  index++;
              else
                while(index<slots.size() and slots[index])
                  index++;
              if(index>=slots.size())
                valid=valid and variadic and not named;
              else if(slots[index])
                valid=false;
              else
                {
                  bool reference=false;
                  if constexpr(requires {parameters[index].isRef;})
                    reference=parameters[index].isRef;
                  slots[index]=reference==argument.reference?argument.type:CompletionBinding{};
                }
            }
      std::vector<CompletionBinding> result;
      for(size_t i=0;i<parameters.size();i++)
        result.push_back(slots[i]?*slots[i]:visit(parameters[i].def));
      
      return valid?std::optional{std::move(result)}:std::nullopt;
    }
    
    /// Infer results from parameter-specialized AST bodies and per-instance struct fields
    CompletionBinding callResult(const CompletionBinding& callable,
                                 const std::vector<Argument>& supplied,
                                 bool& uncertainEffects,
                                 std::set<std::string>& writes)
    {
      CompletionBinding result;
      bool validBody=true;
      if(callable.scriptFunction)
        if(const auto range=script.ranges.find(callable.scriptFunction->body.get());range!=script.ranges.end())
          for(const auto& error : script.diagnostics)
            if(error.range.begin.offset>=range->second.begin.offset and error.range.begin.offset<range->second.end.offset)
              validBody=false;
      if(callable.scriptFunction and validBody and depth<64)
        {
          const auto& parameters=callable.scriptFunction->pars;
          const auto types=arguments(parameters.list,parameters.isVariadic(),supplied);
          CompletionAnalyzer body{.script=script,.remaining=remaining};
          body.returnMode=true;
          body.depth=depth+1;
          body.capturedNames=capturedNames;
          for(const auto& [name,binding] : bindings)
            body.capturedNames.insert(name);
          for(const auto& parameter : parameters.list)
            body.capturedNames.erase(parameter.name);
          for(const auto& [name,binding] : callable.definitionBindings)
            if(const auto current=bindings.find(name);current!=bindings.end() and sameBinding(current->second,*binding))
              body.bindings.emplace(name,current->second);
          if(types)
            for(size_t i=0;i<parameters.list.size();i++)
              body.bindings[parameters.list[i].name]=(*types)[i];
          if(const auto* list=std::get_if<ASTNodesNode>(callable.scriptFunction->body.get()))
            body.statements(*list,false);
          else
            body.visit(callable.scriptFunction->body);
          uncertainEffects=body.invalid;
          writes.insert(body.capturedWrites.begin(),body.capturedWrites.end());
          if(types and not body.invalid and not body.fallsThrough and body.returned)
            result=*body.returned;
        }
      else if(callable.scriptStruct and callable.info.callable and callable.callResult)
        {
          const auto& members=callable.scriptStruct->members.list;
          if(const auto types=arguments(members,false,supplied))
            {
              result=*callable.callResult;
              for(size_t i=0;i<members.size();i++)
                result.memberTypes[members[i].name]=std::make_shared<CompletionBinding>((*types)[i]);
            }
        }
      else if(callable.callResult)
        result=*callable.callResult;
      
      return result;
    }
    
    /// Call signatures and return types come from the named binding or expression callee
    CompletionBinding analyze(const FuncCallNode& node)
    {
      CompletionBinding callable;
      if(matches(node.fun,node.location))
        found=true;
      else if(node.callee)
        callable=visit(node.callee);
      else if(const auto entry=bindings.find(node.fun);entry!=bindings.end())
        callable=entry->second;
      std::vector<Argument> supplied;
      for(const auto& argument : node.args.list)
        {
          if(matches(argument.name,argument.location))
            {
              found=true;
              restricted=true;
              for(const auto& form : callable.info.argumentForms)
                for(auto name : form)
                  {
                    const size_t begin=name.find_first_not_of("[& ");
                    name=begin==std::string::npos?"":name.substr(begin);
                    name=name.substr(0,name.find_first_of(" =]"));
                    if(not name.empty() and name!="..." and name.front()!='<')
                      choices.emplace(name,CompletionBinding{.info={"parameter"},.documentation=callable.documentation});
                  }
            }
          else
            {
              const auto* reference=std::get_if<UnaryRefNode>(argument.expr.get());
              supplied.push_back({argument.name,visit(reference?reference->arg:argument.expr),reference!=nullptr});
            }
        }
      
      bool uncertainEffects=false;
      std::set<std::string> writes;
      const auto result=found?CompletionBinding{}:callResult(callable,supplied,uncertainEffects,writes);
      if(not found)
        {
          for(const auto& name : writes)
            {
              if(bindings.contains(name))
                bindings[name]={};
              if(returnMode and capturedNames.contains(name))
                {
                  invalid=true;
                  capturedWrites.insert(name);
                }
            }
          for(const auto& argument : node.args.list)
            if(const auto* reference=std::get_if<UnaryRefNode>(argument.expr.get()))
              {
                if(const auto* name=std::get_if<IdNode>(reference->arg.get()))
                  {
                    if(bindings.contains(name->name))
                      bindings[name->name]={};
                    if(returnMode and capturedNames.contains(name->name))
                      {
                        capturedWrites.insert(name->name);
                        invalid=true;
                      }
                  }
                else if(std::holds_alternative<MemberAccessNode>(*reference->arg))
                  {
                    for(auto& [name,binding] : bindings)
                      binding.memberTypes.clear();
                    if(returnMode)
                      invalid=true;
                  }
              }
          if(uncertainEffects)
            {
              for(auto& [name,binding] : bindings)
                binding.memberTypes.clear();
              if(returnMode)
                invalid=true;
            }
        }
      
      return result;
    }
    
    /// Retain constructor declarations; specialize member types only at each call
    CompletionBinding analyze(const StructDefNode& node)
    {
      CompletionBinding value;
      value.info={.kind="script struct constructor",.description=node.description,.callable=true,
        .argumentForms=FunctionArgumentForms(1),.location=node.location};
      value.scriptStruct=&node;
      value.callResult=std::make_shared<CompletionBinding>();
      value.callResult->scriptStruct=&node;
      value.callResult->info.kind="script struct instance '"+node.name+"'";
      for(const auto& member : node.members.list)
        {
          value.callResult->members.emplace(member.name,EntityInfo{
            .kind="script member",.description=member.description,.location=member.location});
          value.info.argumentForms[0].push_back(member.def?"["+member.name+" = "+spelling(member.def)+"]":member.name);
          auto child=*this;
          child.visit(member.def);
          if(child.found)
            adopt(child);
        }
      value.documentation=completionDocumentation(node.name,value.info);
      if(not found and node.name!=recoveryMarker and node.name!=completionMarker)
        bindings[node.name]=value;
      
      return value;
    }
    
    /// Add one reachable return, discarding conflicting or unknown result types
    CompletionBinding analyze(const ReturnNode& node)
    {
      const auto value=visit(node.arg);
      if(returnMode)
        {
          mergeReturn(value);
          fallsThrough=false;
        }
      
      return value;
    }
    
    /// Merge one known return into the function's result
    void mergeReturn(const CompletionBinding& value)
    {
      if(not returned)
        returned=value;
      else if(not sameBinding(*returned,value))
        returned=CompletionBinding{};
    }
    
    /// Identify condition effects conservatively using AST alternatives
    bool pureCondition(const std::shared_ptr<ASTNode>& node) const
    {
      bool pure=true;
      if(const auto* binary=std::get_if<BinOpNode>(node.get()))
        pure=pureCondition(binary->arg1) and pureCondition(binary->arg2);
      else if(const auto* unary=std::get_if<UnOpNode>(node.get()))
        pure=pureCondition(unary->arg);
      else
        pure=std::holds_alternative<IdNode>(*node) or std::holds_alternative<ValueNode>(*node);
      
      return pure;
    }
    
    /// Merge types only from paths which continue past the conditional
    CompletionBinding analyze(const IfNode& node)
    {
      visit(node.subNodes[0]);
      if(returnMode and not pureCondition(node.subNodes[0]))
        invalid=true;
      auto left=*this;
      auto right=*this;
      left.returned.reset();
      right.returned.reset();
      left.visit(node.subNodes[1]);
      if(node.subNodes.size()>2)
        right.visit(node.subNodes[2]);
      if(left.found)
        adopt(left);
      else if(right.found)
        adopt(right);
      else if(not found)
        {
          for(auto& [name,value] : bindings)
            {
              if(not returnMode or (left.fallsThrough and right.fallsThrough))
                value=sameBinding(left.bindings.at(name),right.bindings.at(name))?left.bindings.at(name):CompletionBinding{};
              else if(left.fallsThrough)
                value=left.bindings.at(name);
              else if(right.fallsThrough)
                value=right.bindings.at(name);
            }
          if(returnMode)
            {
              invalid=invalid or left.invalid or right.invalid;
              capturedWrites.insert(left.capturedWrites.begin(),left.capturedWrites.end());
              capturedWrites.insert(right.capturedWrites.begin(),right.capturedWrites.end());
              if(left.returned)
                mergeReturn(*left.returned);
              if(right.returned)
                mergeReturn(*right.returned);
              fallsThrough=left.fallsThrough or right.fallsThrough;
            }
        }
      
      return {};
    }
    
    /// Loop-local names remain inside the loop; iteration-dependent types are unknown
    CompletionBinding analyze(const ForNode& node)
    {
      auto body=*this;
      for(const auto& child : node.subNodes)
        body.visit(child);
      if(body.found)
        adopt(body);
      else
        for(auto& [name,value] : bindings)
          if(not sameBinding(value,body.bindings.at(name)))
            value={};
      if(returnMode)
        {
          invalid=true;
          capturedWrites.insert(body.capturedWrites.begin(),body.capturedWrites.end());
        }
      
      return {};
    }
    
    /// Intermediate grammar nodes can survive in a recovered forest
    CompletionBinding analyze(const FuncArgNode& node)
    {
      return visit(node.expr);
    }
    
    /// A parameter default remains an unevaluated AST expression
    CompletionBinding analyze(const FuncParNode& node)
    {
      return visit(node.def);
    }
    
    /// A member default remains an unevaluated AST expression
    CompletionBinding analyze(const StructMemberNode& node)
    {
      return visit(node.def);
    }
    
    /// Visit surviving argument-list fragments after recovery
    CompletionBinding analyze(const FuncArgListNode& node)
    {
      for(const auto& child : node.list)
        analyze(child);
      
      return {};
    }
    CompletionBinding analyze(const FuncParListNode& node)
    {
      for(const auto& child : node.list)
        analyze(child);
      
      return {};
    }
    CompletionBinding analyze(const StructMemberListNode& node)
    {
      for(const auto& child : node.list)
        analyze(child);
      
      return {};
    }
  };

  /// Seed analysis from registered runtime metadata, preserving the original scope
  std::map<std::string,CompletionBinding> environmentBindings(const Environment& environment)
  {
    std::map<std::string,CompletionBinding> result;
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
                binding.scriptStruct=definition->structDefNode;
                binding.callResult=std::make_shared<CompletionBinding>();
                binding.callResult->scriptStruct=definition->structDefNode;
                binding.callResult->members=scriptMembers(*definition->structDefNode);
              }
          }
        result.emplace(name,std::move(binding));
      }
    
    return result;
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
  if(byteOffset<=source.text.size())
    {
      const auto parsed=parseScriptForEditor(source);
      result=completeScriptAt(environment,parsed,byteOffset);
      if(not result.applicable and scriptCompletionRange(parsed,byteOffset))
        result=completeScriptAt(environment,parseScriptForEditor(source,byteOffset),byteOffset);
    }
  
  return result;
}

ScriptCompletion completeScriptAt(const Environment& environment,
                                  const ParsedScript& parsed,
                                  const size_t byteOffset)
{
  ScriptCompletion result;
  const auto& source=parsed.source;
  if(byteOffset<=source.text.size())
    {
      if(const auto target=parsed.completion?parsed.completion:scriptCompletionRange(parsed,byteOffset);
         target and target->begin.offset<=byteOffset and byteOffset<=target->end.offset)
        {
          result.begin=target->begin.offset;
          result.end=target->end.offset;
          CompletionAnalyzer visitor{parsed,environmentBindings(environment),target};
          for(const auto& root : parsed.roots)
            {
              if(const auto range=parsed.ranges.find(root.get());range!=parsed.ranges.end())
                visitor.activeRange=range->second;
              if(const auto* list=std::get_if<ASTNodesNode>(root.get()))
                visitor.statements(*list,false);
              else
                visitor.visit(root);
            }
          result.applicable=visitor.found;
          const auto prefix=std::string_view(source.text).substr(result.begin,byteOffset-result.begin);
          const auto add=[&](const std::string& name,
                            const CompletionBinding& binding)
          {
            if(name.starts_with(prefix) and name!=recoveryMarker and name!=completionMarker)
              result.candidates.push_back({name,binding.info,binding.documentation});
          };
          for(const auto& [name,binding] : visitor.restricted?visitor.choices:visitor.bindings)
            add(name,binding);
          if(not visitor.restricted)
            for(const auto& keyword : parsed.keywords)
              if(not visitor.bindings.contains(keyword))
                add(keyword,CompletionBinding{.info={"keyword"}});
          std::sort(result.candidates.begin(),result.candidates.end(),[](const auto& left,
                                                                       const auto& right)
          {
            return left.name<right.name;
          });
          result.candidates.erase(std::unique(result.candidates.begin(),result.candidates.end(),[](const auto& left,
                                                                                                 const auto& right)
          {
            return left.name==right.name;
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
          writeScriptLispString(output,candidate.name);
          output<<' ';
          writeScriptLispString(output,candidate.info.kind);
          output<<' ';
          writeScriptLispString(output,candidate.documentation);
          output<<' ';
          writeScriptLispString(output,candidate.info.location?candidate.info.location->source:"");
          const SourceCursor cursor=candidate.info.location?candidate.info.location->cursor:SourceCursor{};
          output<<' '<<cursor.line<<' '<<cursor.column<<' '<<cursor.offset<<')';
        }
      output<<"))";
    }
  
  return output.str();
}

std::string formatScriptCompletionForBash(const ScriptCompletion& completion,
                                         const std::string_view text,
                                         const size_t byteOffset)
{
  std::ostringstream output;
  if(completion.applicable)
    {
      output<<completion.begin<<'\0'<<completion.end<<'\0';
      const size_t following=text.find_first_not_of(" \t\r\n",completion.end);
      const bool insertCall=completion.candidates.size()==1 and byteOffset==completion.end and
        (following==std::string_view::npos or text[following]!='(');
      for(const auto& candidate : completion.candidates)
        {
          output<<candidate.name;
          if(insertCall and candidate.info.callable and not candidate.info.argumentForms.empty())
            {
              output<<'(';
              bool first=true;
              for(std::string_view argument : candidate.info.argumentForms.front())
                {
                  // Reference markers describe binding, not an expression to insert
                  if(argument.starts_with('&'))
                    argument.remove_prefix(1);
                  
                  // Only plain names are required: bracketed defaults and ... are omitted
                  const bool identifier=not argument.empty() and
                    ((argument.front()>='a' and argument.front()<='z') or
                     (argument.front()>='A' and argument.front()<='Z') or argument.front()=='_') and
                    std::all_of(argument.begin(),argument.end(),[](const char c)
                    {
                      return (c>='a' and c<='z') or (c>='A' and c<='Z') or
                        (c>='0' and c<='9') or c=='_';
                    });
                  if(identifier)
                    {
                      output<<(first?"":", ")<<argument;
                      first=false;
                    }
                }
              output<<')';
            }
          output<<'\0';
        }
    }
  
  return output.str();
}
