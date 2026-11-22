#ifndef _NASTY_HPP
#define _NASTY_HPP

/**
  General-purpose script grammar, abstract syntax tree, and evaluator for NISSA.
  Applications supply domain-specific host types, functions, and execution.

  getParser() recognizes source text. getParseTreeExecuctor() maps parser
  actions to ASTNode values, which Evaluator executes through std::visit.
  Value stores runtime data; Environment gives names shared storage cells.
  HostFunction exposes C++ calls, while HObject and HObjectInfo expose concrete
  C++ values through registered member accessors and printing functions.

  Keep the AST alive while using script functions or struct definitions:
  their runtime representations retain pointers into the tree.
*/

#include <any>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <optional>
#include <set>
#include <tuple>
#include <utility>
#include <type_traits>
#include <typeindex>
#include <limits>

#include <parsePact.hpp>
#include <tuples/tupleHasType.hpp>

#include "hostMembers.hpp"
#include "scriptSource.hpp"

using namespace pp::internal;

/// Build the ParsePact parser for the script grammar and its semantic actions
pp::Grammar getParser();

/// Extract immediately preceding /// and //! lines from a token's skipped text
std::string precedingDocumentation(const pp::internal::ShiftNode& token);

/// Assignment operations shared by script values and host member setters
/// The enumerator order matches the integer tags emitted by the grammar
enum class AssignMode{PLAIN,PROD,DIV,SUM,SUB};

struct Evaluator;

struct Environment;

struct ASTNodesNode;

struct ValueNode;

struct PrePostfixOpNode;

struct FuncNode;

struct FuncDefNode;

struct FuncCallNode;

struct ForNode;

struct IfNode;

struct UnOpNode;

struct UnaryRefNode;

struct BinOpNode;

struct IdNode;

struct AssignNode;

struct ReturnNode;

struct FuncParNode;

struct FuncArgNode;

struct FuncArgListNode;

struct FuncParListNode;

struct StructMemberNode;

struct StructMemberListNode;

struct StructDefNode;

struct SubscribeNode;

struct MemberAccessNode;

/// All executable and intermediate syntax nodes produced by parser actions
/// Child nodes use shared pointers to allow recursive variant definitions
using ASTNode=
  std::variant<ASTNodesNode,
	       UnaryRefNode,
	       UnOpNode,
	       BinOpNode,
	       ForNode,
	       IfNode,
	       IdNode,
	       PrePostfixOpNode,
	       FuncNode,
	       FuncDefNode,
	       FuncCallNode,
	       FuncParNode,
	       FuncParListNode,
	       SubscribeNode,
	       FuncArgNode,
	       FuncArgListNode,
	       StructMemberNode,
	       StructMemberListNode,
	       StructDefNode,
	       MemberAccessNode,
	       ReturnNode,
	       ValueNode,
	       AssignNode>;

struct ValueRef;

struct Function;

struct HostFunction;

struct HObjectInfo;

struct ValuesList;

struct StructDef;

struct Struct;

struct HObject;

struct FuncArg;

struct Value;

struct HObjectMemberRef;

/// Reference to a script storage element
struct ValueRef
{
  std::shared_ptr<Value> ref;
};

/// Sequence of shared value cells, used for script lists and variadic arguments
/// Copying the list copies the cell handles, not the contained values
struct ValuesList
{
  std::vector<std::shared_ptr<Value>> data;
  
  Value& operator[](const int& i);
};

/// Evaluated call argument; an empty name denotes a positional argument
struct FuncArg
{
  std::string name;
  
  std::shared_ptr<Value> value;
};

/// Script function together with the environment captured at definition time
/// funcNode borrows the AST; env keeps the captured scope alive
struct Function
{
  const FuncNode* funcNode;
  
  std::shared_ptr<Environment> env;
  
  /// Original declaration name, or empty for a lambda
  std::string name{};
};

/// Alternative argument lists for display; brackets denote optional arguments
using FunctionArgumentForms=
  std::vector<std::vector<std::string>>;

/// Static result description; an absent C++ type means a dynamically determined value
struct ValueTypeInfo
{
  /// Concrete scalar or host payload type, when known
  std::optional<std::type_index> cppType{};
  
  /// Registered host description, weak to avoid cycles through binary overloads
  std::weak_ptr<const HObjectInfo> hostType{};
};

/// C++ callable exposed to scripts, receiving prepared arguments and the evaluator
/// An empty argumentNames signature permits positional arguments only
struct HostFunction :
  std::function<Value(std::vector<std::shared_ptr<Value>>&,Evaluator&)>
{
  /// Argument names used for binding; constructor arguments may be omitted
  std::vector<std::string> argumentNames;
  
  /// Whether this callable constructs a registered host type
  bool isConstructor{false};
  
  /// Optional explanation for verbose discovery
  std::string description{};
  
  /// Display-only raw signatures; otherwise derive named arguments or report an unknown signature
  FunctionArgumentForms argumentForms{};
  
  /// Declared result type, available to tooling without invoking this callable
  ValueTypeInfo resultType{};
};

/// Callable script struct definition borrowing its declaration from the AST
struct StructDef
{
  const StructDefNode* structDefNode;
};

/// Script struct instance with member cells ordered as in its AST declaration
/// Copying an instance shares its member cells; structDefNode borrows the AST
struct Struct
{
  const StructDefNode* structDefNode;
  
  std::vector<std::shared_ptr<Value>> members;
  
  inline Struct(const StructDef& def,
		std::vector<std::shared_ptr<Value>>& members);
};

inline std::ostream& operator<<(std::ostream& os,
				const Struct& s);

struct HObjectInfo;

/// Runtime alternatives; monostate represents an uninitialized or absent value
/// Integers use int or int64_t when larger; boolean results remain int
using ValueT=
  std::variant<std::monostate,std::string,int,int64_t,double,Function,FuncArg,HostFunction,ValueRef,ValuesList,StructDef,Struct,HObject>;

/// Concrete C++ value stored in std::any, together with its shared type metadata
/// Copying copies the value according to its C++ copy semantics; metadata is shared
struct HObject
{
  std::shared_ptr<const HObjectInfo> info;
  
  std::any data;
  
  inline HObjectMemberRef getMemberRef(const std::string& memberName);
};

inline std::ostream& operator<<(std::ostream& os,
				const HObject& obj);

/// Runtime value wrapper providing variant assignment and optional conversion
/// Copying follows the semantics of the active alternative, including shared cells
struct Value
{
  /// Returns the stored value as T, or an empty optional if no type matches
  /// First tries to retrieve T itself from the variant. Otherwise, tries the
  /// proposed From types in order and converts the first matching alternative
  /// to T using static_cast. Each From must be a variant alternative and
  /// support that cast
  /// T need not be a variant alternative; HObject values are handled separately
  template <typename T,
	    VariantAlternative<ValueT>...From>
  std::optional<T> tryAs() const
  {
    // Holds the converted value only when a proposed source type matches
    std::optional<T> result;
    
    // Try to retrieve the stored value as S and convert it to T
    auto attempt=
      [&]<typename S>()
      {
	return std::visit([&](const auto& source)
	{
	  // Match the proposed source type, not any type convertible to T
	  if constexpr(std::is_same_v<S,std::remove_cvref_t<decltype(source)>>)
	    {
	      result.emplace(static_cast<T>(source));
	      return true;
	    }
	  else
	    return false;
	},data);
      };
    
    // Try T first, then From in order, stopping at the first match
    (attempt.template operator()<T>() or ... or
     attempt.template operator()<From>());
    
    return result;
  }
  
  Value& operator=(const ValueT& oth)
  {
    data=oth;
    
    return *this;
  }
  
  // ValueT operator*(const Value& oth) const
  // {
  //   return {std::visit([](const auto& lhs,
  // 			  const auto& rhs)->ValueT
  //   {
  //      if constexpr(requires {lhs * rhs;})
  // 	 return {lhs*rhs};
  //      else
  // 	 errorEmitter("Cannot "  "assign the types: ",typeid(lhs).name()," and ",typeid(rhs).name());
       
  //      return {};
  //   },this->data,oth.data)};
  // }
  
  ValueT data;
};

/// Recover the actual C++ type T from the given positional host-function argument
/// Named script arguments have already been placed in their signature positions
/// T is either a ValueT alternative or the concrete payload of an HObject
/// Require an exact type match, except int64_t also accepts int by value
/// Other returned references borrow the argument cells and must not outlive them
template <typename T>
decltype(auto) hostArgumentTypeConcretizer(const std::vector<std::shared_ptr<Value>>& arguments,
                                    const std::size_t index,
                                    const char* functionName)
{
  if(index>=arguments.size() or not arguments[index])
    errorEmitter("missing argument ",index," to ",functionName);
  
  if constexpr(std::is_same_v<T,int64_t>)
    {
      const auto value=arguments[index]->tryAs<int64_t,int>();
      
      if(not value)
        errorEmitter("invalid argument ",index," to ",functionName," (expected integer)");
      
      return int64_t(*value);
    }
  else
    {
      const ValueT& data=
        arguments[index]->data;
      const T* value=nullptr;
      
      if constexpr(VariantAlternative<T,ValueT>)
        value=std::get_if<T>(&data);
      else if(const auto* object=std::get_if<HObject>(&data))
        value=std::any_cast<T>(&object->data);
      
      if(not value)
        errorEmitter("invalid argument ",index," to ",functionName,
                     " (expected ",typeid(T).name(),")");
      
      return *value;
    }
}

/// Recover positional arguments as a tuple of their actual C++ types for std::apply
/// Types follow signature order; every argument is required and extras are rejected
/// Tuple elements borrow argument cells, except widened integers which are owned
template <typename...T>
auto hostArgumentsTypeConcretizer(const std::vector<std::shared_ptr<Value>>& arguments,
                                                   const char* functionName)
{
  if(arguments.size()>sizeof...(T))
    errorEmitter("too many arguments to ",functionName);
  
  return [&]<std::size_t...I>(std::index_sequence<I...>)
  {
    return std::tuple<decltype(hostArgumentTypeConcretizer<T>(arguments,I,functionName))...>{
      hostArgumentTypeConcretizer<T>(arguments,I,functionName)...};
  }(std::index_sequence_for<T...>{});
}

inline std::ostream& operator<<(std::ostream& os,
				const Value& v)
{
  std::visit([&os](const auto& v)
  {
    if constexpr(requires {os<<v;})
      os<<v;
    else
      os<<"unprintable type: "<<typeid(v).name()<<"\n";
  },v.data);
  
  return os;
}

/// Read/write handle for a host member, implemented through paired functions
/// Unlike ValueRef, it need not refer to an actual Value cell
struct HObjectMemberRef
{
  using Getter=
    std::function<Value()>;
  
  using Setter=
    std::function<void(const AssignMode&,const Value&)>;
  
  using Accessors=
    std::pair<Getter,Setter>;
  
  Accessors accessors;
  
  Value get() const
  {
    return accessors.first();
  }
  
  void set(const AssignMode& mode,
	   const Value& value)
  {
    accessors.second(mode,value);
  }
};

/// Enable runtime arithmetic checks in script operators and native compound assignments
inline constexpr bool performArithmeticValidation=true;

/// Reject invalid integer quotients before division or remainder is evaluated
/// Check the promoted operand type, preserving mixed and floating-point arithmetic
template <typename L,
          typename R>
void validateIntegerDivision(const L& lhs,
                             const R& rhs)
{
  if constexpr(performArithmeticValidation and
               std::is_integral_v<L> and std::is_integral_v<R>)
    {
      using Result=
        decltype(lhs/rhs);
      const Result left=lhs;
      const Result right=rhs;
      
      if(right==0)
        errorEmitter("integer divisor is zero");
      else if constexpr(std::is_signed_v<Result>)
        if(left==std::numeric_limits<Result>::min() and right==-1)
          errorEmitter("integer division overflow");
    }
}

/// Adapt a native binary functor and a diagnostic name to script dispatch
/// std::array is a structural C++20 template argument for the name
template <typename NativeOperator,
          auto Name>
struct BinaryOperator
{
  static constexpr const char* name=Name.data();
  
  template <typename L,
            typename R>
  requires requires(const L& lhs,
                    const R& rhs)
  {
    NativeOperator{}(lhs,rhs);
  }
  auto operator()(const L& lhs,
                  const R& rhs) const
  {
    if constexpr(std::is_same_v<NativeOperator,std::divides<>> or
                 std::is_same_v<NativeOperator,std::modulus<>>)
      validateIntegerDivision(lhs,rhs);
    
    return NativeOperator{}(lhs,rhs);
  }
};

/// Addition with its script error name
using SumOperator=
  BinaryOperator<std::plus<>,std::to_array("sum")>;

/// Subtraction with its script error name
using DifferenceOperator=
  BinaryOperator<std::minus<>,std::to_array("difference")>;
  
/// Multiplication with its script error name
using ProductOperator=
  BinaryOperator<std::multiplies<>,std::to_array("product")>;
  
/// Division with its script error name
using DivisionOperator=
  BinaryOperator<std::divides<>,std::to_array("division")>;

/// Integer remainder with its script error name
using RemainderOperator=
  BinaryOperator<std::modulus<>,std::to_array("remainder")>;

/// Shared description of a registered host type and its supported operations
/// Member providers bind accessors to an individual value supplied at access time
struct HObjectInfo
{
  /// Script-visible type name, also used in diagnostics
  std::string name;
  
  /// Concrete payload identity used when selecting registered binary overloads
  std::type_index cppType{typeid(void)};
  
  /// Member accessors and their static value type, sharing one registration entry
  struct MemberProvider :
    std::function<HObjectMemberRef::Accessors(std::any&)>
  {
    /// Callable interface retained when replacing or instrumenting a provider
    using Function=
      std::function<HObjectMemberRef::Accessors(std::any&)>;
    
    /// Replacing the callable preserves its registered value type
    using Function::operator=;
    
    /// Script-visible member type; scalar conversions match the getter
    ValueTypeInfo valueType{};
  };
  
  /// Factories binding named member accessors to the supplied concrete value
  std::unordered_map<std::string,MemberProvider> memberAccessorsProvider;
  
  /// Callable and static result description for one concrete pair of binary operands
  struct BinaryFunction :
    std::function<Value(const std::any&,const std::any&)>
  {
    /// Result metadata retained beside the callable, without invoking it
    ValueTypeInfo resultType{};
  };
  
  /// Concrete overloads for one binary operator
  struct BinaryOverloads
  {
    /// Right operand types when this host object is on the left
    std::unordered_map<std::type_index,BinaryFunction> rightOperands;
    
    /// Left scalar types when this host object is on the right
    std::unordered_map<std::type_index,BinaryFunction> leftOperands;
  };
  
  /// Binary operators keyed by their operation type
  std::unordered_map<std::type_index,BinaryOverloads> binaryOperations;
  
  /// Write the value to the stream
  std::function<void(std::ostream&,const std::any&)> streamOp;
};

/// Dispatch a registered host operation or apply it to script values directly
/// Host-host and host-scalar pairs preserve their original operand order
template <typename Operator>
Value applyBinaryOperator(const Value& lhs,
                          const Value& rhs)
{
  return std::visit([]<typename L,
                     typename R>(const L& left,
                                 const R& right)->Value
  {
    const auto typeName=[](const auto& operand)->std::string
    {
      if constexpr(std::is_same_v<std::remove_cvref_t<decltype(operand)>,HObject>)
        return operand.info->name;
      else
        return typeid(operand).name();
    };
    
    if constexpr(std::is_same_v<L,HObject>)
      {
        const auto operation=left.info->binaryOperations.find(typeid(Operator));
        if(operation!=left.info->binaryOperations.end())
          {
            if constexpr(std::is_same_v<R,HObject>)
              {
                const auto found=operation->second.rightOperands.find(right.data.type());
                if(found!=operation->second.rightOperands.end())
                  return found->second(left.data,right.data);
              }
            else
              {
                const auto found=operation->second.rightOperands.find(typeid(R));
                if(found!=operation->second.rightOperands.end())
                  return found->second(left.data,std::any(right));
              }
          }
      }
    else if constexpr(std::is_same_v<R,HObject>)
      {
        const auto operation=right.info->binaryOperations.find(typeid(Operator));
        if(operation!=right.info->binaryOperations.end())
          {
            const auto found=operation->second.leftOperands.find(typeid(L));
            if(found!=operation->second.leftOperands.end())
              return found->second(std::any(left),right.data);
          }
      }
    else if constexpr(requires {Operator{}(left,right);})
      return {Operator{}(left,right)};
    
    errorEmitter("unsupported ",Operator::name,": ",typeName(left),
                 " and ",typeName(right));
    return {};
  },lhs.data,rhs.data);
}

/// Bind a native binary operation to host metadata
/// A scalar operand has null metadata; result metadata must already exist
template <typename Operator,
          typename L,
          typename R>
requires requires(const L& lhs,
                  const R& rhs)
{
  Operator{}(lhs,rhs);
}
void registerHObjectBinaryOperator(const std::shared_ptr<HObjectInfo>& lhsInfo,
                                   const std::shared_ptr<HObjectInfo>& rhsInfo,
                                   const std::shared_ptr<HObjectInfo>& resultInfo)
{
  if(not resultInfo or (not lhsInfo and not rhsInfo))
    errorEmitter("a host binary operation requires operand and result metadata");
  
  HObjectInfo::BinaryFunction function={
    [result=std::weak_ptr<HObjectInfo>(resultInfo)](const std::any& lhs,
                                                    const std::any& rhs)->Value
  {
    auto info=result.lock();
    if(not info)
      errorEmitter("binary operation result metadata is no longer registered");
    
    return {HObject{info,Operator{}(std::any_cast<const L&>(lhs),
                                    std::any_cast<const R&>(rhs))}};
  },{resultInfo->cppType,resultInfo}};
  
  if(lhsInfo)
    lhsInfo->binaryOperations[typeid(Operator)].rightOperands[typeid(R)]=
      std::move(function);
  else
    rhsInfo->binaryOperations[typeid(Operator)].leftOperands[typeid(L)]=
      std::move(function);
}

inline std::ostream& operator<<(std::ostream& os,
				const HObject& obj)
{
  if(const auto& sOp=obj.info->streamOp)
    sOp(os,obj.data);
  else
    errorEmitter("Unable to print host object of type \"",obj.info->name,"\" (stream operator not registered)");
  
  return os;
}


/// Bind accessors to this value, reporting an unknown member name as an error
inline HObjectMemberRef HObject::getMemberRef(const std::string& memberName)
{
  auto& ap=
    info->memberAccessorsProvider;
  
  auto rt=
    ap.find(memberName);
  
  if(rt==ap.end())
    errorEmitter("Trying to subscribe member \"",memberName,"\" but the subscribed object does not have such member");
  
  return {rt->second(this->data)};
}

/////////////////////////////////////////////////////////////////

/// Variable name, resolved against the evaluator's current environment
struct IdNode
{
  std::string name;
  
  /// Documentation used when this identifier defines a variable
  std::string description{};
  
  /// Position of the identifier token when parsed from source
  std::optional<ScriptLocation> location{};
};

/// Call argument expression with an optional name for named argument binding
struct FuncArgNode
{
  std::string name;
  
  std::shared_ptr<ASTNode> expr;
};

/// Call arguments in source order, before evaluation and parameter binding
struct FuncArgListNode
{
  std::vector<FuncArgNode> list;
};

/// List indexing expression: base[subscr]
struct SubscribeNode
{
  std::shared_ptr<ASTNode> base;
  
  std::shared_ptr<ASTNode> subscr;
};

/// Named member access on a script struct or registered host object
struct MemberAccessNode
{
  std::shared_ptr<ASTNode> base;
  
  std::string member;
};

/// Literal or parser token represented as a runtime Value
struct ValueNode
{
  Value value;
  
  /// Documentation before a parser token, consumed only by definitions
  std::string description{};
  
  /// Original token position, owned independently of the parse tree
  std::optional<ScriptLocation> location{};
  
  /// Explicitly initialize the nested runtime value without aggregate brace elision
  template <typename V>
  static ValueNode fromValue(const V& value)
  {
    return {.value=Value{value}};
  }
  
  /// Convert a shifted token while discarding ordinary skipped comments
  static ValueNode fromShift(const pp::internal::ShiftNode& token)
  {
    const auto position=token.location.value_or(pp::internal::SourceLocation{{},token.cursor});
    
    return {.value=Value{std::string(token.txt)},.description=precedingDocumentation(token),
            .location=ScriptLocation{position.source.empty()?"<script>":std::string(position.source),position.cursor}};
  }
};

/// Statement block, evaluated in a child environment in source order
struct ASTNodesNode
{
  std::vector<std::shared_ptr<ASTNode>> list;
};

/// Conditional with children ordered as condition, then branch, optional else branch
struct IfNode
{
  std::vector<std::shared_ptr<ASTNode>> subNodes;
};

/// Function parameter with a name, reference flag, and optional default expression
struct FuncParNode
{
  std::string name;
  
  bool isRef;
  
  std::shared_ptr<ASTNode> def;
  
  /// Position of a named parameter's declaration
  std::optional<ScriptLocation> location{};
  
  static constexpr const char VARIADIC_PAR_NAME[]="__VA__ARGS__";
  
  /// Build a variadic parameter whose reference policy applies to each element
  static FuncParNode getVariadic(const bool isRef)
  {
    return {.name=VARIADIC_PAR_NAME,.isRef=isRef,.def=std::make_shared<ASTNode>(ValueNode{.value{ValuesList{}}})};
  }
};

/// Feature-provider to D, to find the passed name in the "list" member of D
template <typename D>
struct Searchable
{
  /// In case of missing, list.size() is returned
  auto iParOfName(const std::string& name) const
  {
    const auto& list=
      ((const D*)this)->list;
    
    auto it=
      std::find_if(list.begin(),
		   list.end(),
		   [&name](const auto& par)
		   {
		     return par.name==name;
		   });
    
    return std::distance(list.begin(),it);
  }
};

/// Ordered function parameters, optionally ending in the reserved variadic parameter
struct FuncParListNode:
  Searchable<FuncParListNode>
{
  std::vector<FuncParNode> list;
  
  bool isVariadic() const
  {
    return list.size() and list.back().name==FuncParNode::VARIADIC_PAR_NAME;
  }
  
  /// Append a variadic parameter with an empty list as its default
  void makeVariadic(const bool isRef)
  {
    if(isVariadic())
      errorEmitter("Cannot make a list of pars variadic twice");
    
    list.push_back(FuncParNode::getVariadic(isRef));
    
  };
  
  /// Build a signature consisting only of a variadic parameter
  static FuncParListNode getEmptyVariadic(const bool isRef)
  {
    FuncParListNode list;
    
    list.makeVariadic(isRef);
    
    return list;
  }
};

/// Loop with children ordered as initializer, condition, increment, and body
struct ForNode
{
  std::vector<std::shared_ptr<ASTNode>> subNodes;
};

/// Function parameters and body; evaluating this node creates a closure
struct FuncNode
{
  FuncParListNode pars;
  
  std::shared_ptr<ASTNode> body;
  
  /// Immediately preceding documentation, retained by callable aliases through this AST
  std::string description{};
  
  /// Position of fun or lambda, retained by aliases
  std::optional<ScriptLocation> location{};
};

/// Named function declaration, installing a closure in the current environment
struct FuncDefNode
{
  std::string name;
  
  std::shared_ptr<FuncNode> fun;
};

/// Call to a named script function, host function, or script struct constructor
struct FuncCallNode
{
  std::string fun;
  
  FuncArgListNode args;
  
  /// Position of the called name
  std::optional<ScriptLocation> location{};
};

/// Declared script struct member with an optional default expression
struct StructMemberNode
{
  std::string name;
  
  std::shared_ptr<ASTNode> def;
  
  /// Position of the member's declaration
  std::optional<ScriptLocation> location{};
};

/// Ordered script struct members with lookup by member name
struct StructMemberListNode:
  Searchable<StructMemberListNode>
{
  std::vector<StructMemberNode> list;
};

/// Named script struct declaration; member order also defines positional arguments
struct StructDefNode
{
  std::string name;
  
  StructMemberListNode members;
  
  /// Immediately preceding documentation shared by constructor aliases
  std::string description{};
  
  /// Position of the struct keyword
  std::optional<ScriptLocation> location{};
};

inline std::ostream& operator<<(std::ostream& os,
				const Struct& s)
{
  const StructDefNode& def=*s.structDefNode;
  
  os<<"Struct "<<def.name<<"\n";
  for(size_t i=0;i<def.members.list.size();i++)
    os<<" "<<def.members.list[i].name<<" "<<*s.members[i]<<"\n";
  
  return os;
};

inline Struct::Struct(const StructDef& def,
		      std::vector<std::shared_ptr<Value>>& members) :
  structDefNode(def.structDefNode),
  members(members)
{
  const std::string& name=
    def.structDefNode->name;
  
  diagnostic("Creating struct ",name,"\n");
    
  const size_t nArgs=members.size();
  const size_t nMembers=
    def.structDefNode->members.list.size();
  
  if(members.size()!=nMembers)
    errorEmitter("Calling struct ",name," constructor with ",nArgs," args when expecting ",nMembers);
}

/// Plain or compound assignment to a writable expression
struct AssignNode
{
  std::shared_ptr<ASTNode> lhs;
  
  std::shared_ptr<ASTNode> rhs;
  
  AssignMode mode;
};

/// Increment/decrement expression using the passed operator
struct PrePostfixOpNode
{
  std::shared_ptr<ASTNode> arg;
  
  std::function<Value(Value&)> op;
};

/// Explicit reference expression used to pass a script value cell by reference
struct UnaryRefNode
{
  std::shared_ptr<ASTNode> arg;
};

/// Unary expression using the passed operator
struct UnOpNode
{
  std::shared_ptr<ASTNode> arg;
  
  std::function<Value(const Value&)> op;
};

/// Binary expression with the passed operator
struct BinOpNode
{
  /// Whether both operands or only the logically necessary operands are evaluated
  enum class EvaluationMode
  {
    /// Evaluate both operands and invoke the callable
    EAGER,
    
    /// Skip the right operand when the left operand is false
    LOGICAL_AND,
    
    /// Skip the right operand when the left operand is true
    LOGICAL_OR
  };
  
  /// Left operand, evaluated before the right operand
  std::shared_ptr<ASTNode> arg1;
  
  /// Right operand, evaluated only when required by the mode
  std::shared_ptr<ASTNode> arg2;
  
  /// Callable used by eager operations; logical modes produce Boolean results directly
  std::function<Value(const Value&,const Value&)> op;
  
  /// Operand evaluation policy selected when constructing the expression
  EvaluationMode evaluationMode{EvaluationMode::EAGER};
};

/// Returning value - node
struct ReturnNode
{
  std::shared_ptr<ASTNode> arg;
};

/// Discovery metadata shared by classification and verbose listings
struct EntityInfo
{
  /// Human-readable entity kind
  std::string kind;
  
  /// Documentation supplied by its native registration
  std::string description{};
  
  /// Whether this value can be called from the script
  bool callable{false};
  
  /// Argument forms for callable values, in declaration order
  FunctionArgumentForms argumentForms{};
  
  /// Script declaration position, absent for native entities
  std::optional<ScriptLocation> location{};
};

/// Describe native signatures without changing positional or named argument binding
FunctionArgumentForms entityArgumentForms(const HostFunction& function);

/// Describe script parameters, including references, defaults, and variadic arguments
FunctionArgumentForms entityArgumentForms(const Function& function);

/// Describe required and defaulted members of a script constructor
FunctionArgumentForms entityArgumentForms(const StructDef& definition);

/// Inspect every runtime value alternative; new alternatives require an overload
struct EntityInfoVisitor
{
  /// Reject unhandled types, including types convertible to an existing alternative
  template <typename T>
  EntityInfo operator()(const T&) const=delete;
  
  // Define one exact visitor overload; the deleted fallback still rejects new types
#define PROVIDE_ENTITY_INFO(TYPE,NAME,...) \
  EntityInfo operator()(const TYPE& NAME) const \
  { \
    return {__VA_ARGS__}; \
  }
  
  /// Uninitialized values
  PROVIDE_ENTITY_INFO(std::monostate,,"uninitialized variable")
  
  /// Strings
  PROVIDE_ENTITY_INFO(std::string,,"variable (string)")
  
  /// Integers
  PROVIDE_ENTITY_INFO(int,,"variable (int)")
  
  /// Wide integers
  PROVIDE_ENTITY_INFO(int64_t,,"variable (int64)")
  
  /// Real numbers
  PROVIDE_ENTITY_INFO(double,,"variable (double)")
  
  /// Script functions and anonymous lambdas
  PROVIDE_ENTITY_INFO(Function,function,function.name.empty()?"script lambda":"embedded function '"+function.name+"'",function.funcNode->description,true,entityArgumentForms(function),function.funcNode->location)
  
  /// Named argument values
  PROVIDE_ENTITY_INFO(FuncArg,,"variable (named argument)")
  
  /// Native functions and constructors
  PROVIDE_ENTITY_INFO(HostFunction,function,function.isConstructor?"host struct constructor":"host function",function.description,true,entityArgumentForms(function))
  
  /// References to variable cells
  PROVIDE_ENTITY_INFO(ValueRef,,"variable (reference)")
  
  /// Lists of values
  PROVIDE_ENTITY_INFO(ValuesList,,"variable (list)")
  
  /// Script struct constructors
  PROVIDE_ENTITY_INFO(StructDef,definition,"script struct constructor",definition.structDefNode->description,true,entityArgumentForms(definition),definition.structDefNode->location)
  
  /// Script struct instances
  PROVIDE_ENTITY_INFO(Struct,object,"script struct instance '"+object.structDefNode->name+"'")
  
  /// Native object instances
  PROVIDE_ENTITY_INFO(HObject,object,"host object '"+object.info->name+"'")
  
#undef PROVIDE_ENTITY_INFO
};

/// Select the visible bindings included in a discovery listing
enum class EntitySelection{Functions,Variables,All};

/// Format visible bindings; verbose mode groups variables and separates documented signatures
std::string formatEntityListing(const Environment& environment,
                                const EntitySelection selection,
                                const bool verbose);

/// Format the language overview or documentation for one visible entity
/// An empty keyword selects the overview; unknown names produce a diagnostic
std::string formatHelp(const Environment& environment,
                       const std::string& keyword="");

/// Contains shared variable cells and possibly the enclosing scope
/// Creates only if no larger scope has the name
struct Environment
{
  std::shared_ptr<Environment> parent;
  
  std::unordered_map<std::string,std::shared_ptr<Value>> varTable;

  /// Descriptions attached to bindings in this scope
  std::unordered_map<std::string,std::string> descriptions;
  
  /// First script declaration position of each binding in this scope
  std::unordered_map<std::string,ScriptLocation> locations;
  
  /// List selected visible names alphabetically, respecting inner bindings
  std::vector<std::string> entityNames(const EntitySelection selection=EntitySelection::All) const
  {
    std::set<std::string> seen;
    std::vector<std::string> names;
    for(const Environment* scope=this;scope!=nullptr;scope=scope->parent.get())
      for(const auto& [name,value] : scope->varTable)
        if(seen.insert(name).second)
          {
            const bool callable=value and std::visit(EntityInfoVisitor{},value->data).callable;
            
            if(selection==EntitySelection::All or
               (selection==EntitySelection::Functions and callable) or
               (selection==EntitySelection::Variables and not callable))
              names.push_back(name);
          }
    
    std::sort(names.begin(),names.end());
    
    return names;
  }
  
  /// List the visible callable names
  std::vector<std::string> functionNames() const
  {
    return entityNames(EntitySelection::Functions);
  }
  
  /// Describe a visible binding without changing its value
  void describe(const std::string& name,
                const std::string& description)
  {
    Environment* scope=this;
    while(scope!=nullptr and not scope->varTable.contains(name))
      scope=scope->parent.get();
    
    if(scope==nullptr)
      errorEmitter("Cannot describe unknown entity ",name);
    
    scope->descriptions[name]=description;
  }
  
  /// Inspect a binding in its nearest scope, applying binding-specific documentation
  EntityInfo entityInfo(const std::string& name) const
  {
    const Environment* scope=this;
    while(scope!=nullptr and not scope->varTable.contains(name))
      scope=scope->parent.get();
    
    if(scope==nullptr)
      errorEmitter("Cannot inspect unknown entity ",name);
    
    const auto& value=scope->varTable.at(name);
    EntityInfo info=
      value?std::visit(EntityInfoVisitor{},value->data):EntityInfo{"uninitialized variable"};
    
    if(const auto found=scope->descriptions.find(name);found!=scope->descriptions.end())
      info.description=found->second;
    
    if(not info.location)
      if(const auto found=scope->locations.find(name);found!=scope->locations.end())
        info.location=found->second;
    
    return info;
  }
  
  /// Explain a binding's kind and optional description in its nearest scope
  std::string entityDetails(const std::string& name) const
  {
    const EntityInfo info=entityInfo(name);
    std::string details=info.kind;
    
    if(not info.description.empty())
      details+=" - "+info.description;
    
    return details;
  }
  
  /// Find the outermost environment in the parent chain
  Environment& uppermost()
  {
    Environment* res=this;
    
    while(Environment* par=&*res->parent)
      res=par;
    
    return *res;
  }
  
  /// Look up a name in this scope and its parents; return null when absent
  std::shared_ptr<Value> find(const std::string& name)
  {
    if(auto it=varTable.find(name);it!=varTable.end())
      return it->second;
    else
      if(parent)
	return parent->find(name);
      else
	return {};
  }
  
  /// Requires an existing value, reports an error if the name cannot be found
  Value& at(const std::string& name)
  {
    std::shared_ptr<Value> f=find(name);
    
    if(f==nullptr)
      errorEmitter("Trying to use uninitialized variable ",name);
    
    return *f;
  }
  
  /// Returns existing variable or insert an uninitialized cell, returning true if the cell was not existing
  std::pair<std::shared_ptr<Value>,bool> getRefOrInsert(const std::string& name)
  {
    if(auto f=find(name))
      return {f,false};
    else
      return {varTable[name]=std::make_shared<Value>(std::monostate{}),true};
  }
  
  /// Access a visible binding, creating an uninitialized local cell when needed
  Value& operator[](const std::string& name)
  {
    return *getRefOrInsert(name).first;
  }
  
  // void print(const size_t& i=0)
  // {
  //   using namespace pp::internal;
  
  //   diagnostic("----- ",i," -----\n");
  //   for(const auto& [name,v] : varTable)
  //     std::visit([&name](const auto& v)
  //     {
  // 	if constexpr(Streamable<decltype(v)>)
  // 	  diagnostic(name,"=",v,"\n");
  // 	else
  // 	  diagnostic(name," of unprintable type: ",typeid(decltype(v)).name(),"\n");
  //     },*v);
    
  //   if(parent)
  //     parent->print(i+1);
  // }
  
  Environment(std::shared_ptr<Environment> parent=nullptr) :
    parent(parent)
  {
  }
};

/// Function activation environment whose parent is the captured definition scope
struct CallFrame
{
  std::shared_ptr<Environment> env;
  
  CallFrame(std::shared_ptr<Environment> parentEnv=nullptr) :
    env{std::make_shared<Environment>(parentEnv)}
  {
  }
};

/// Extract the expected AST alternative for a parser action, diagnosing a mismatch
/// comm supplies optional context for the diagnostic; subNode must be non-null
template <typename T>
T& fetch(std::shared_ptr<ASTNode>& subNode,
	 const char* comm=nullptr)
{
  T* s=
    std::get_if<T>(&*subNode);
  
  if(not s)
    errorEmitter("subNode: ",(comm?comm:"")," is not of the required type ",typeid(T).name()," but is of type ",variantInnerTypeName(*subNode));
  
  return *s;
}

/// Extract a parser action child by index and require the expected AST alternative
template <typename T>
T& fetch(std::vector<std::shared_ptr<ASTNode>>& list,
	 const size_t& i)
{
  using namespace pp::internal;
  
  if(const size_t n=list.size();n<i)
    errorEmitter(n," nodes received, aksed for node #",i);
  
  return fetch<T>(list[i],std::to_string(i).c_str());
}

/// Start a list node of type L with an element of AST alternative T
template <typename T,
	  typename L>
std::shared_ptr<ASTNode> makeFirstOfList(std::shared_ptr<ASTNode>& elem)
{
  return std::make_shared<ASTNode>(L{.list{fetch<T>(elem)}});
}

/// Move an AST element into an existing list node and return that list handle
template <typename T,
	  typename L>
std::shared_ptr<ASTNode>& appendToList(std::shared_ptr<ASTNode>& list,
				       std::shared_ptr<ASTNode>& elem)
{
 fetch<L>(list).list.push_back(std::move(fetch<T>(elem)));
 
 return list;
}

/// Append the second list's elements by moving them into the first list
template <typename L>
std::shared_ptr<ASTNode>& mergeLists(std::shared_ptr<ASTNode>& list1,
				     std::shared_ptr<ASTNode>& list2)
{
  auto& first=fetch<L>(list1).list;
  auto& second=fetch<L>(list2).list;
  
  first.insert(first.end(),std::make_move_iterator(second.begin()),std::make_move_iterator(second.end()));
  
  return list1;
}

/// Associate grammar action names with AST construction functions
pp::ParseTreeExecutor<ASTNode,ValueNode> getParseTreeExecuctor(const std::vector<std::string_view>& requiredActions);

/// Apply the requested C++ assignment expression when it is well formed
/// Compound assignments retain native C++ promotion and conversion rules
/// Unsupported operand combinations report an interpreter error
template <typename L,
	  typename R>
void trySet(L& lhs,
	    const AssignMode& mode,
	    const R& rhs)
{
  if(mode==AssignMode::DIV)
    validateIntegerDivision(lhs,rhs);
  
 switch(mode)
   {
#define CASE(MODE,SYMBOL,DESCR)						\
     case AssignMode::MODE:						\
       if constexpr(requires {lhs SYMBOL rhs;})				\
	 lhs SYMBOL rhs;						\
       else								\
	 errorEmitter("Cannot " DESCR "assign the types: ",typeid(lhs).name()," and ",typeid(rhs).name()); \
       break
     
      CASE(PLAIN,=,"");
      CASE(PROD,*=,"multiply-");
      CASE(DIV,/=,"divide-");
      CASE(SUM,+=,"sum-");
      CASE(SUB,-=,"sub-");
      
#undef CASE
    }
}

/// Assign script values while preserving native compound-assignment conversions
/// Plain assignment may replace the stored type; host arithmetic uses registration
inline void trySet(Value& lhs,
                   const AssignMode& mode,
                   const Value& rhs)
{
  if(mode==AssignMode::PLAIN)
    lhs=rhs;
  else if(std::holds_alternative<HObject>(lhs.data) or
          std::holds_alternative<HObject>(rhs.data))
    switch(mode)
      {
#define CASE(MODE,OPERATOR) \
      case AssignMode::MODE: \
        lhs=applyBinaryOperator<OPERATOR>(lhs,rhs); \
        break
        
      CASE(SUM,SumOperator);
      CASE(SUB,DifferenceOperator);
      CASE(PROD,ProductOperator);
      CASE(DIV,DivisionOperator);
      
#undef CASE
      
      default:
	errorEmitter("unsupported compound assignment for host objects");
      }
  else
    std::visit([&mode](auto& left,
		       const auto& right)
    {
      trySet(left,mode,right);
    },lhs.data,rhs.data);
}

/// Register the built-in functions and constants in a new evaluator
void registerEvaluatorBuiltins(Evaluator& evaluator);

/// AST visitor with scope management, call binding, and runtime value access
struct Evaluator
{
  /// Result of evaluating an expression: a reference, a temporary, or an empty result
  /// Statement-only visitors may return an empty result; consume values only from expressions
  struct EvalResult
  {
    using Ref=
      std::variant<std::monostate,ValueRef,HObjectMemberRef>;
    
    Ref ref;
    
    /// Whether this result carries a script cell or host member reference
    bool hasRef() const
    {
      return not std::holds_alternative<std::monostate>(ref);
    }
    
    /// Temporary value, absent when evaluation produces no value
    std::optional<Value> value;
    
    /// Whether a temporary value is available independently of any reference
    bool hasValue() const
    {
      return value.has_value();
    }
    
    /// Whether execution must unwind to the enclosing script function
    bool isReturn{false};
    
    /// Assign through either reference kind; reject results that are not writable
    void set(const AssignMode& mode,
	     const Value& v)
    {
      std::visit(pp::internal::Overload{[&v,
			   &mode](ValueRef& vr)
      {
	if(mode!=AssignMode::PLAIN)
	  ensureInitialized(*vr.ref);
	
	trySet(*vr.ref,mode,v);
      },
	    [&v,
	     &mode](HObjectMemberRef& hr)
	    {
	      hr.set(mode,v);
	    },
	    [](const std::monostate&)
	    {
	      errorEmitter("trying to set a non lhs");
	    }},ref);
    }
    
    /// Require a concrete script Value cell, for operations needing direct mutation
    /// Host member accessors do not expose such a cell and must use set() instead
    Value& asLhs()
    {
      ValueRef* _cRefP=std::get_if<ValueRef>(&ref);
      
      if(not _cRefP)
	errorEmitter("trying to use a non-concrete type as lhs");
      
      ValueRef& cRefP=*_cRefP;
      
      if(not cRefP.ref)
	errorEmitter("trying to access an evaluation result as an lhs when it is not");
      
      return *cRefP.ref;
    }
    
    /// Reject a Value containing the uninitialized monostate alternative
    static void ensureInitialized(const Value& v)
    {
      if(std::get_if<std::monostate>(&v.data))
	errorEmitter("using uninitialized value");
    }
    
    /// Read the current referent, or return the stored temporary value
    /// Reference reads are deferred so later writes are visible when this is called
    Value asRhs() const
    {
      if(hasRef())
	return
	  std::visit(pp::internal::Overload{[this](const ValueRef& vr)->Value
	  {
	    if(const std::shared_ptr<Value>& r=vr.ref)
	      {
		ensureInitialized(*r);
		
		return *r;
	      }
	    else
	      {
		if(not hasValue())
		  errorEmitter("trying to read an empty evaluator result");
		
		ensureInitialized(*value);
		
		return *value;
	      }
	  },
		[](const HObjectMemberRef& hr)->Value
		{
		  return hr.get();
		},
		[](const std::monostate&) ->Value
		{
		  errorEmitter("impossible");
		  
		  return {std::monostate{}};
		}},ref);
      else
	{
	  if(not hasValue())
	    errorEmitter("trying to read an empty evaluator result");
	  
	  return *value;
	}
    }
  };
  
  /// Convert a condition or logical operand using the script's Boolean conversion rules
  static bool asBool(const Value& value)
  {
    return std::visit([](const auto& data)
    {
      if constexpr(std::is_convertible_v<decltype(data),bool>)
        return static_cast<bool>(data);
      else
        errorEmitter("Cannot convert the type to bool");
      
      return false;
    },value.data);
  }
  
  /// Retained activation frames, starting with the root environment
  std::vector<CallFrame> callFrames;
  
  /// Retain the latest metadata per C++ type, including types used only as results
  std::unordered_map<std::type_index,std::shared_ptr<HObjectInfo>> hostObjectTypes;
  
  /// Scope used for name lookup and insertion during the current evaluation
  std::shared_ptr<Environment> curEnv;
  
  /// Temporarily change the current scope and restore it on destruction
  struct EnvironmentWindUnwinder
  {
    Evaluator& ev;
    
    std::shared_ptr<Environment> prevEnv;
    
    EnvironmentWindUnwinder(Evaluator& ev,
			    std::shared_ptr<Environment> newEnv) :
      ev{ev},prevEnv{ev.curEnv}
    {
      ev.curEnv=newEnv;
    }
    
    ~EnvironmentWindUnwinder()
    {
      ev.curEnv=prevEnv;
    }
  };
  
  /// Enter a child scope; retain the returned guard for the lifetime of that scope
  EnvironmentWindUnwinder pushNewEnvironment()
  {
    return {*this,std::make_shared<Environment>(curEnv)};
  }
  
  /// Create the root scope and install constants, list utilities, printing, and math calls
  Evaluator()
  {
    curEnv=callFrames.emplace_back().env;
    registerEvaluatorBuiltins(*this);
  }
  
  /// Construct a script struct from an argument vector already ordered by member
  EvalResult createStruct(const StructDef& def,
			  std::vector<std::shared_ptr<Value>>& args)
  {
    const std::string& name=def.structDefNode->name;
    
    diagnostic("Creating struct ",name,"\n");
    
    const size_t nArgs=args.size();
    const size_t nMembers=def.structDefNode->members.list.size();
    
    if(args.size()!=nMembers)
      errorEmitter("Calling struct ",name," constructor with ",nArgs," args when expecting ",nMembers);
      
    Struct res(def,args);
    
    return {.value=Value{res}};
  }
  
  /// Execute a script closure with arguments already arranged in parameter order
  EvalResult callEmbeddedFunction(const Function& f,
				  std::vector<std::shared_ptr<Value>>& args)
  {
    const FuncNode& fn=*f.funcNode;
    const FuncParListNode& pars=fn.pars;
    
    const size_t nArgs=args.size();
    const size_t nPars=pars.list.size();
    
    if(nArgs!=nPars)
      errorEmitter("Calling ",(pars.isVariadic()?"variadic ":" "),"function with ",nArgs," args when expecting ",nArgs);
    
    CallFrame& cf=callFrames.emplace_back(f.env);
    
    EnvironmentWindUnwinder windUnwind(*this,cf.env);
    
    for(size_t iPar=0;iPar<nPars;iPar++)
      {
	const FuncParNode& par=pars.list[iPar];
	
	std::shared_ptr<Value> arg=args[iPar];
	
	if(par.isRef)
	  {
	    const auto referent=[&par](const std::shared_ptr<Value>& argument)
	    {
	      const auto* reference=
		std::get_if<ValueRef>(&argument->data);
	      
	      if(not reference or not reference->ref)
		errorEmitter("parameter ",par.name," has to be taken by ref but passed argument is not a ref");
	      
	      return reference->ref;
	    };
	    
	    if(pars.isVariadic() and iPar+1==nPars)
	      {
		const auto* arguments=
		  std::get_if<ValuesList>(&arg->data);
		
		if(not arguments)
		  errorEmitter("variadic parameter must receive an argument list");
		
		// Bind the pack's elements to caller cells without changing the supplied list
		ValuesList references;
		references.data.reserve(arguments->data.size());
		for(const auto& argument : arguments->data)
		  references.data.push_back(referent(argument));
		
		arg=std::make_shared<Value>(std::move(references));
	      }
	    else
	      arg=referent(arg);
	  }
	
	if(not curEnv->varTable.try_emplace(par.name,arg).second)
	  errorEmitter("argument ",par.name," already specified");
        
        if(par.location)
          curEnv->locations.emplace(par.name,*par.location);
      }
    
    EvalResult result=
      std::visit(*this,*fn.body);
    
    // Consume any explicit return signal at the function-call boundary
    result.isReturn=false;
    
    return result;
  }
  
  /// Dispatch prepared arguments to a script closure, host function, or struct constructor
  EvalResult callFunction(const Value& vFun,
			  std::vector<std::shared_ptr<Value>>& args)
  {
    if(const Function* f=std::get_if<Function>(&vFun.data))
      return callEmbeddedFunction(*f,args);
    else
      if(const HostFunction* hf=std::get_if<HostFunction>(&vFun.data))
	return {.value=(*hf)(args,*this)};
      else
	if(const StructDef* sd=std::get_if<StructDef>(&vFun.data))
	  return createStruct(*sd,args);
	else
	  errorEmitter("Variable is of type ",variantInnerTypeName(vFun.data)," not supported by the callFunction backend");
    
    return {};
  }
  
  /// Find a callable by name and invoke it with already prepared arguments
  EvalResult callFunction(const std::string& name,
			  std::vector<std::shared_ptr<Value>>& args)
  {
    diagnostic("Going to call function ",name.c_str(),"\n");
    
    return callFunction(curEnv->at(name),args);
  }
  
  /// Return the literal as a temporary value
  EvalResult operator()(const ValueNode& valueNode)
  {
    return {.value=valueNode.value};
  }
  
  /// Resolve or create a variable cell and return its writable reference
  EvalResult operator()(const IdNode& idNode)
  {
    auto [cell,inserted]=curEnv->getRefOrInsert(idNode.name);
    if(inserted and idNode.location)
      curEnv->locations.emplace(idNode.name,*idNode.location);
    
    return {.ref=ValueRef(cell)};
  }
  
  /// Evaluate an explicit reference expression using a concrete script cell
  EvalResult operator()(const UnaryRefNode& uRNode)
  {
    EvalResult tmp=std::visit(*this,*uRNode.arg);
    
    ValueRef* _vr=std::get_if<ValueRef>(&tmp.ref);
    if(not _vr)
      errorEmitter("cannot take the reference of a non-embedded object");
    
    ValueRef& vr=*_vr;
    if(not vr.ref)
      errorEmitter("subNode has not returned a reference");
      
    return {.value=Value{vr}};
  }
  
  /// Install a callable struct definition, rejecting a name already visible in scope
  EvalResult operator()(const StructDefNode& structDefNode)
  {
    const std::string& name=
      structDefNode.name;
    
    if(curEnv->find(name))
      errorEmitter("Redefining a struct which has a name ",name," already defined");
    
    (*curEnv)[name]=StructDef{.structDefNode=&structDefNode};
    
    return {};
  }
  
  /// Index a list with a nonnegative integer, checking its bounds
  EvalResult operator()(const SubscribeNode& subscribeNode)
  {
    EvalResult _b=std::visit(*this,*subscribeNode.base);
    const Value s=std::visit(*this,*subscribeNode.subscr).asRhs();
    Value& b=_b.hasRef()?_b.asLhs():*_b.value;
    
    if(ValuesList* vl=std::get_if<ValuesList>(&b.data))
      if(const auto mi=s.tryAs<int64_t,int>())
	{
	  const int64_t i=*mi;
	  
	  if(i<0)
	    errorEmitter("trying to subscribe with negative index ",i);
	  
	  if(i>=vl->data.size())
	    errorEmitter("trying to subscribe withindex ",i," greater than list size");
	  
	  std::shared_ptr<Value> r=vl->data[i];
	  
	  if(_b.hasRef())
	    return {.ref{ValueRef{r}}};
	  else
	    return {.value=*r};
	}
      else
	errorEmitter("index of subscrition is not an integer but is of type: ",variantInnerTypeName(s.data));
    else
      errorEmitter("trying to subscribe non-list of type: ",variantInnerTypeName(b.data));
    
    return {};
  }
  
  /// Resolve a script or host member, preserving writability when the base has a reference
  EvalResult operator()(const MemberAccessNode& memberAccessNode)
  {
    EvalResult _b=std::visit(*this,*memberAccessNode.base);
    Value b=_b.asRhs();
    
    const std::string& memberName=memberAccessNode.member;
    
    if(Struct* es=std::get_if<Struct>(&b.data))
      {
	const StructDefNode* sd=es->structDefNode;
	
	size_t iMember=sd->members.iParOfName(memberName);
	
	if(iMember==es->members.size())
	  errorEmitter("trying to access to member ",memberName," not present in embedded struct of type: ",sd->name);
	
	std::shared_ptr<Value> mf=es->members[iMember];
	if(_b.hasRef())
	  return {.ref={ValueRef{mf}}};
	else
	  return {.value=*mf};
      }
    else
      if(HObject* hObject=std::get_if<HObject>(&b.data))
	{
          // Validate the member even when this expression is only used as an lhs
          HObjectMemberRef member=hObject->getMemberRef(memberName);
          if(not _b.hasRef())
            return {.value=member.get()};
	  
          // Retain the parent reference, not a pointer into a temporary value
          // Read it again on each access so RHS changes to siblings are preserved
          HObjectMemberRef nested{{
            [base=_b,memberName]()->Value
            {
              Value parent=base.asRhs();
              return unvariant<HObject>(parent.data).getMemberRef(memberName).get();
            },
            [base=_b,memberName](const AssignMode& mode,const Value& value) mutable
            {
              Value parent=base.asRhs();
              unvariant<HObject>(parent.data).getMemberRef(memberName).set(mode,value);
              base.set(AssignMode::PLAIN,parent);
            }}};
          return {.ref{nested}};
        }
      else
	errorEmitter("trying to access member ",memberName," of a non-stuct of type: ",variantInnerTypeName(b.data));
    
    return {};
  }
  
  /// Execute the loop in a child scope, checking its condition each iteration
  EvalResult operator()(const ForNode& forNode)
  {
    EnvironmentWindUnwinder windUnwind=pushNewEnvironment();
    
    for(std::visit(*this,*forNode.subNodes[0]);
        asBool(std::visit(*this,*forNode.subNodes[1]).asRhs());
	std::visit(*this,*forNode.subNodes[2]))
      {
        EvalResult result=
          std::visit(*this,*forNode.subNodes[3]);
        
        // Leave the loop before its increment and forward the return to the function
        if(result.isReturn)
          return result;
      }
    
    return {};
  }
  
  /// Evaluate the condition in a child scope and execute the selected branch
  EvalResult operator()(const IfNode& ifNode)
  {
    EnvironmentWindUnwinder windUnwind=
      pushNewEnvironment();
    
    int iNodeToRun=0;
    
    if(asBool(std::visit(*this,*ifNode.subNodes[0]).asRhs()))
      iNodeToRun=1;
    else
      if(ifNode.subNodes.size()>=3)
	iNodeToRun=2;
    
    if(iNodeToRun)
      return std::visit(*this,*ifNode.subNodes[iNodeToRun]);
    else
      return {};
  }
  
  /// Parameter declarations are consumed during call binding, not executed directly
  EvalResult operator()(const FuncParNode&)
  {
    return {};
  }
  
  /// Parameter lists are metadata; direct evaluation produces no value
  EvalResult operator()(const FuncParListNode&)
  {
    return {};
  }
  
  /// Evaluate a member default, or produce monostate when no default is present
  EvalResult operator()(const StructMemberNode& structMemberNode)
  {
    if(std::shared_ptr<ASTNode> def=structMemberNode.def)
      return {.value=Value{std::visit(*this,*def).asRhs().data}};
    else
      return {.value=Value{std::monostate{}}};
  }
  
  /// Evaluate member defaults in declaration order into a list of value cells
  EvalResult operator()(const StructMemberListNode& structMemberListNode)
  {
    ValuesList memberList;
    
    const size_t nMembers=structMemberListNode.list.size();
    
    memberList.data.reserve(nMembers);
    
    diagnostic("Preparing a list of ",nMembers," members\n");
    
    for(const StructMemberNode& structMemberNode : structMemberListNode.list)
      memberList.data.push_back(std::make_shared<Value>((*this)(structMemberNode).asRhs()));
    
    return {.value=Value{memberList}};
  }
  
  /// Evaluate an argument expression and retain its optional argument name
  EvalResult operator()(const FuncArgNode& argNode)
  {
    FuncArg arg{.name=argNode.name,
		.value=std::make_shared<Value>(std::visit(*this,*argNode.expr).asRhs())};

    return {.value=Value{std::move(arg)}};
  }
  
  /// Evaluate call arguments in source order before matching them to parameters
  EvalResult operator()(const FuncArgListNode& argListNode)
  {
    ValuesList argList;
    
    const size_t nArgs=argListNode.list.size();
    
    argList.data.reserve(nArgs);
    
    diagnostic("Preparing a list of ",nArgs," args\n");
    
    for(const FuncArgNode& argNode : argListNode.list)
      argList.data.push_back(std::make_shared<Value>((*this)(argNode).asRhs()));
    
    return {.value=Value{argList}};
  }
  
  /// Prepare named and positional arguments according to the callable's signature
  /// Script defaults are evaluated in the caller's scope. Named host signatures
  /// leave omitted slots null so the host constructor can retain its C++ defaults
  EvalResult operator()(const FuncCallNode& funcCallNode)
  {
    Value __args=(*this)(funcCallNode.args).asRhs();
    
    ValuesList* _args=std::get_if<ValuesList>(&__args.data);
    if(not _args)
      errorEmitter("FuncCallNode args did not evaluate to a ValuesList but of type ",variantInnerTypeName(__args.data));
    
    const std::string fName=
      funcCallNode.fun;
    
    diagnostic("Going to call function ",fName,"\n");
    
    std::vector<std::shared_ptr<Value>> args;
    if(not curEnv->find(fName) and funcCallNode.location)
      funcCallNode.location->error("Trying to call unknown function ",fName);
    
    const Value& _fun=curEnv->at(fName);
    
    /// Fill a parameter slot once, using the same diagnostic for every callable
    const auto insertArgument=
      [&args,
       &fName](const size_t index,
               const std::string& name,
               const std::shared_ptr<Value>& value)
      {
        if(args[index])
          errorEmitter("argument ",name," specified twice to ",fName);
        else
          args[index]=value;
      };
    
    if([[maybe_unused]] const HostFunction* fun=std::get_if<HostFunction>(&_fun.data))
      {
        const auto& names=
	  fun->argumentNames;
        if(not names.empty())
          args.resize(names.size());
	
        size_t positional{};
	for(std::shared_ptr<Value> _arg : _args->data)
	  {
	    FuncArg& arg=unvariant<FuncArg>(_arg->data);
	    
            if(names.empty())
              {
                if(arg.name!="")
                  errorEmitter("trying to call a host function with named arg ",arg.name);
                args.push_back(arg.value);
              }
            else
              {
                size_t i=positional;
                if(arg.name.empty())
                  positional++;
                else
                  {
                    i=0;
                    while(i<names.size() and names[i]!=arg.name)
                      i++;
                    if(i==names.size())
                      errorEmitter("unknown argument ",arg.name," to ",fName);
                  }
                if(i>=args.size())
                  errorEmitter("too many arguments to ",fName);
                insertArgument(i,names[i],arg.value);
              }
	  }
      }
    else
      if(const Function* fun=std::get_if<Function>(&_fun.data))
	{
	  const FuncNode& fn=*fun->funcNode;
	  const FuncParListNode& pars=fn.pars;
	  
	  const size_t nPars=pars.list.size();
	  const bool isVariadic=pars.isVariadic();
	  
	  args.resize(nPars);
	  if(isVariadic)
	    args.back()=std::make_shared<Value>(ValuesList{});
	  
	  auto insertNamedOrUnnamed=
	    [isVariadic,
	     &args,
	     &_args,
	     &pars,
	     &fName,
	     &insertArgument](const bool insertNamed)
	    {
	      for(const std::shared_ptr<Value>& __arg : _args->data)
		{
		  FuncArg& _arg=unvariant<FuncArg>(__arg->data);
		  std::shared_ptr<Value> arg=_arg.value;
		  
		  const std::string& aName=_arg.name;
		  
		  const bool isNamedArg=aName!="";
		  
		  if(not (isNamedArg xor insertNamed))
		    {
		      size_t iPar{};
		      
		      if(isNamedArg)
			{
			  iPar=
			    pars.iParOfName(aName);
			  
			  diagnostic("Plugging named argument: ",aName," as arg n.",iPar," of function ",fName,"\n");
			  
			  if(iPar==pars.list.size())
			    errorEmitter("trying to pass argument ",aName," not expected by the function ",fName);
			}
		      else
			while(iPar<args.size() and args[iPar])
			  iPar++;
		      
		      if(iPar>=args.size())
			{
			  if(isVariadic)
			    std::get_if<ValuesList>(&args.back()->data)->data.push_back(arg);
			  else
			    errorEmitter("passing too many arguments to non-variadic function ",fName);
			}
		      else
			insertArgument(iPar,pars.list[iPar].name,arg);
		    }
		}
	    };
	  
	  for(const int& namedUnnamed : {0,1})
	    insertNamedOrUnnamed(namedUnnamed);
	  
	  for(size_t iPar=0;iPar<pars.list.size();iPar++)
	    if(not args[iPar])
	      {
		const FuncParNode& par=pars.list[iPar];
		
		if(std::shared_ptr<ASTNode> optDef=par.def)
		  args[iPar]=std::make_shared<Value>(std::visit(*this,*optDef).asRhs());
		else
		  errorEmitter("parameter \"",par.name,"\" with no default value unspecified when calling the function");
	      }
	}
      else
	if(const StructDef* _sd=std::get_if<StructDef>(&_fun.data))
	  {
	    const StructDefNode& sd=*_sd->structDefNode;
	    const StructMemberListNode& members=sd.members;
	    
	    const size_t nMembers=members.list.size();
	    args.resize(nMembers);
	    
	    auto insertNamedOrUnnamed=
	      [&args,
	       &_args,
	       &members,
	       &fName,
	       &insertArgument](const bool insertNamed)
	      {
		for(const std::shared_ptr<Value>& __arg : _args->data)
		  {
		    FuncArg& _arg=unvariant<FuncArg>(__arg->data);
		    std::shared_ptr<Value> arg=_arg.value;
		    
		    const std::string& aName=_arg.name;
		    
		    const bool isNamedArg=aName!="";
		    
		    if(not (isNamedArg xor insertNamed))
		      {
			size_t iMember{};
			
			if(isNamedArg)
			  {
			    iMember=
			      members.iParOfName(aName);
			    
			    diagnostic("Plugging named argument: ",aName," as member n.",iMember," of struct ",fName,"\n");
			    
			    if(iMember==members.list.size())
			      errorEmitter("trying to pass argument ",aName," not expected by the struct ",fName);
			  }
			else
			  while(iMember<args.size() and args[iMember])
			    iMember++;
			
			if(iMember>=args.size())
			  errorEmitter("passing too many arguments to struct ",fName);
			else
			  insertArgument(iMember,members.list[iMember].name,arg);
		      }
		  }
	      };
	    
	    for(const int& namedUnnamed : {0,1})
	      insertNamedOrUnnamed(namedUnnamed);
	    
	    for(size_t iMember=0;iMember<members.list.size();iMember++)
	      if(not args[iMember])
		{
		  const StructMemberNode& member=members.list[iMember];
		  
		  if(std::shared_ptr<ASTNode> optDef=member.def)
		    {
		      diagnostic("Member ",member.name," missing value, adding its default\n");
		      
		      args[iMember]=std::make_shared<Value>(std::visit(*this,*optDef).asRhs());
		    }
		  else
		    errorEmitter("member \"",member.name,"\" with no default value unspecified when calling the struct constructor");
		}
	  }
	else
	  errorEmitter("trying to call a non-callable variable of kind: ",variantInnerTypeName(_fun.data));
    
    return callFunction(_fun,args);
  }
  
  /// Bind a named closure, rejecting a name already visible in the environment
  EvalResult operator()(const FuncDefNode& funcDefNode)
  {
    const std::string& name=
      funcDefNode.name;
    
    if(curEnv->find(name))
      errorEmitter("Redefining a function which has a name ",name," already defined");
    
    Value function=(*this)(*funcDefNode.fun).asRhs();
    std::get<Function>(function.data).name=name;
    (*curEnv)[name]=std::move(function);
    
    return {};
  }
  
  /// Capture the current environment alongside a borrowed pointer to the function AST
  EvalResult operator()(const FuncNode& funcNode)
  {
    return
      {.value=Value{Function{.funcNode=&funcNode,
		       .env=curEnv}}};
  }
  
  /// Execute a block in a child scope until a statement signals a return
  EvalResult operator()(const ASTNodesNode& astNodesNode)
  {
    EnvironmentWindUnwinder windUnwind=pushNewEnvironment();
    
    for(const auto& n : astNodesNode.list)
      {
        EvalResult result=
          std::visit(*this,*n);
        
        // Stop this block and forward the return through its enclosing statements
        if(result.isReturn)
          return result;
      }
    
    return {};
  }
  
  /// Evaluate the expression and signal a return through enclosing statements
  EvalResult operator()(const ReturnNode& returnNode)
  {
    EvalResult result=
      std::visit(*this,*returnNode.arg);
    
    // Mark the explicit return so enclosing statements stop executing
    result.isReturn=true;
    
    return result;
  }
  
  /// Apply increment/decrement to a concrete script cell and return the operation result
  EvalResult operator()(const PrePostfixOpNode& prePostfixOpNode)
  {
    return{.value=prePostfixOpNode.op(std::visit(*this,*prePostfixOpNode.arg).asLhs())};
  }
  
  /// Resolve the left-hand reference, evaluate the RHS, then assign through the reference
  EvalResult operator()(const AssignNode& assignNode)
  {
    EvalResult _lhs=std::visit(*this,*assignNode.lhs);
    
    _lhs.set(assignNode.mode,std::visit(*this,*assignNode.rhs).asRhs());
    
    if(assignNode.mode==AssignMode::PLAIN)
      if(const auto* id=std::get_if<IdNode>(assignNode.lhs.get());id and not id->description.empty())
        curEnv->describe(id->name,id->description);
    
    return _lhs;
  }
  
  /// Evaluate the operand and invoke the unary operation function
  EvalResult operator()(const UnOpNode& unOpNode)
  {
    return {.value=unOpNode.op(std::visit(*this,*unOpNode.arg).asRhs())};
  }
  
  /// Evaluate left first and skip the right operand when the logical mode permits
  EvalResult operator()(const BinOpNode& binOpNode)
  {
    const Value left=
      std::visit(*this,*binOpNode.arg1).asRhs();
    
    if(binOpNode.evaluationMode!=BinOpNode::EvaluationMode::EAGER)
      {
        const bool condition=
          asBool(left);
        
        if((binOpNode.evaluationMode==BinOpNode::EvaluationMode::LOGICAL_AND and not condition) or
           (binOpNode.evaluationMode==BinOpNode::EvaluationMode::LOGICAL_OR and condition))
          return {.value=Value{condition}};
        
        return {.value=Value{asBool(std::visit(*this,*binOpNode.arg2).asRhs())}};
      }
    
    const Value right=
      std::visit(*this,*binOpNode.arg2).asRhs();
    
    return {.value=binOpNode.op(left,right)};
  }
};

/// Register metadata without creating a script constructor, for example for results
/// The evaluator keeps this metadata available while scripts are evaluated
template <typename T>
std::shared_ptr<HObjectInfo> registerHObjectType(Evaluator& ev,
						 const std::string& name)
{
  auto info=
    std::make_shared<HObjectInfo>();
  
  info->name=name;
  info->cppType=typeid(T);
  ev.hostObjectTypes.insert_or_assign(typeid(T),info);
  
  return info;
}

/// Registers a value-owning constructor and its fields in the evaluator
/// Tuple order determines positional arguments; omitted fields retain T{} defaults
/// Metadata is shared, but each HObject contains its own copy of T
/// Returns the metadata needed to register T as a member of another host type
/// T must support T{} construction and copying into std::any
/// int, int64_t, double, bool, size_t, and std::string members use scalar accessors; other members
/// require their own metadata. bool and size_t reads become script integers;
/// size_t values outside the int range report an error on read
/// The optional description documents the constructor; otherwise a generic one is supplied
/// The caller may attach additional functions, such as streamOp, to the result
template <typename T,
	  typename...Members>
std::shared_ptr<HObjectInfo> registerHObject(Evaluator& ev,
					     const std::string& name,
					     const std::tuple<Members...>& members,
                                             const std::string& description={})
{
  /// Metadata is retained by the evaluator, including for later product dispatch
  auto info=registerHObjectType<T>(ev,name);
  
  // Constructor argument names, in the same order as the member tuple
  std::vector<std::string> names;
  
  // Register each member and the functions used to access it
  std::apply([&](const auto&...field)
  {
    auto add=
      [&](auto field)
      {
	using M=
	  std::remove_cvref_t<decltype(std::declval<T>().*(field.pointer))>;
	
	if(field.name.empty() or info->memberAccessorsProvider.contains(field.name))
	  errorEmitter("empty or duplicate host member name in ",name,": ",field.name);
	
	/// Types which supports scalar
	using ScalarMemberTypes=
	  std::tuple<int,int64_t,double,bool,size_t,std::string>;
	
	/// Check whether the type is in the list
	constexpr bool scalar=
	  nissa::tupleHasType<ScalarMemberTypes,M>;
	
	// Reuse previously registered metadata for nested host members
	if constexpr(not scalar)
	  if(not field.memberType)
	    {
	      const auto nested=ev.hostObjectTypes.find(typeid(M));
	      if(nested==ev.hostObjectTypes.end())
		errorEmitter("missing host metadata for member ",field.name," in ",name);
	      field.memberType=nested->second;
	    }
	
	names.push_back(field.name);
	
	// Retain exactly the type returned by the getter, including integer normalization
	using ScriptType=
	  std::conditional_t<std::is_same_v<M,bool> or std::is_same_v<M,size_t>,int,M>;
	const ValueTypeInfo memberType{typeid(ScriptType),scalar?nullptr:field.memberType};
	
	// Bind accessors to the value of the particular object being accessed
	info->memberAccessorsProvider.emplace(field.name,
					      HObjectInfo::MemberProvider{[field](std::any& data)
					      {
						// Read a scalar value or wrap a copy of the nested host member
						auto getter=[&data,field]()->Value
						{
						    const M& value=
						      std::any_cast<const T&>(data).*(field.pointer);
						    Value result;
						    
						    if constexpr(scalar)
						      {
							// Script integers also represent booleans and size_t index fields
							if constexpr(std::is_same_v<M,size_t>)
							  if(value>static_cast<size_t>(std::numeric_limits<int>::max()))
							    errorEmitter("host member ",field.name," exceeds the script integer range");
							result=Value{ScriptType(value)};
						      }
						    else
						      result=Value{HObject{field.memberType,value}};
						    
						    return result;
						};
						
						// Apply plain or compound assignment to the concrete C++ member
						auto setter=[&data,field](const AssignMode& mode,
									  const Value& value)
						{
						  M& target=
						    std::any_cast<T&>(data).*(field.pointer);
						  
						  // Host objects assign their payload when its type matches the member
						  if(const HObject* object=std::get_if<HObject>(&value.data))
						    {
						      if(const M* payload=std::any_cast<M>(&object->data))
							trySet(target,mode,*payload);
						      else
							errorEmitter("invalid value for host member ",field.name,
								     " (expected ",typeid(M).name(),")");
						    }
						  else if constexpr(scalar)
						    {
						      // Preserve ordinary script scalar conversion rules
						      std::visit([&](const auto& rhs)
						      {
							trySet(target,mode,rhs);
						      },value.data);
						    }
						  else
						    errorEmitter("invalid value for host member ",field.name,
								 " (expected ",typeid(M).name(),")");
						};
						
						return HObjectMemberRef::Accessors{getter,setter};
					      },memberType});
      };
    
    (add(field),...);
  },members);
  
  // Register the script-visible constructor under the requested type name
  // The evaluator uses names to place positional and named arguments in args
  (*ev.curEnv)[name]=HostFunction{
    [info,names,name](std::vector<std::shared_ptr<Value>>& args,
		      Evaluator&)->Value
    {
      if(args.size()>names.size())
	errorEmitter("too many arguments to ",name);
      
      // Start with T{} defaults; omitted arguments have null entries
      HObject object{info,T{}};
      
      // Initialize supplied members through the same setters as assignment
      for(size_t i=0;i<args.size();i++)
	if(args[i])
	  object.getMemberRef(names[i]).set(AssignMode::PLAIN,*args[i]);
      
      return {object};
    },names,true,description.empty()?"Construct "+name+" from its named members":description,{},{typeid(T),info}};
  
  return info;
}

/// Extract a typed lambda's argument types and invoke it with script values
template <typename Signature>
struct HostCallableSignature;

template <typename Owner,
          typename Return,
          typename...Arguments>
struct HostCallableSignature<Return (Owner::*)(Arguments...) const>
{
  static constexpr std::size_t arity=sizeof...(Arguments);
  
  /// Native return type before conversion to a script value
  using Result=
    std::remove_cvref_t<Return>;
  
  /// Concretize required arguments in signature order, then call the lambda
  template <typename Callable>
  static Value invoke(Callable& callable,
                      const std::vector<std::shared_ptr<Value>>& values,
                      const char* functionName,
                      Evaluator& evaluator)
  {
    auto arguments=
      hostArgumentsTypeConcretizer<std::remove_cvref_t<Arguments>...>(values,functionName);
    
    if constexpr(std::is_void_v<Return>)
      {
        std::apply(callable,arguments);
        
        return {};
      }
    else if constexpr(std::is_same_v<Result,Value>)
      return std::apply(callable,arguments);
    else if constexpr(std::is_constructible_v<ValueT,Return>)
      return {std::apply(callable,arguments)};
    else
      {
        const auto found=evaluator.hostObjectTypes.find(typeid(Result));
        if(found==evaluator.hostObjectTypes.end())
          errorEmitter("unregistered host result type in ",functionName);
        
        return {HObject{found->second,std::apply(callable,arguments)}};
      }
  }
};

/// Inject the calling evaluator without exposing it as a script argument
template <typename Owner,
          typename Return,
          typename...Arguments>
struct HostCallableSignature<Return (Owner::*)(Evaluator&,Arguments...) const> :
  HostCallableSignature<Return (Owner::*)(Arguments...) const>
{
  /// Bind the evaluator, then reuse ordinary script argument concretization
  template <typename Callable>
  static Value invoke(Callable& callable,
                      const std::vector<std::shared_ptr<Value>>& values,
                      const char* functionName,
                      Evaluator& evaluator)
  {
    auto bound=[&](Arguments...arguments)->Return
    {
      return callable(evaluator,std::forward<Arguments>(arguments)...);
    };
    
    return HostCallableSignature<Return (Owner::*)(Arguments...) const>::invoke(
      bound,values,functionName,evaluator);
  }
};

/// Mutable lambdas have the same script argument contract
template <typename Owner,
          typename Return,
          typename...Arguments>
struct HostCallableSignature<Return (Owner::*)(Arguments...)> :
  HostCallableSignature<Return (Owner::*)(Arguments...) const>
{
};

/// Describe a native result without constructing it or evaluating a callback
/// Erased runtime wrappers deliberately leave the result unknown
/// Non-scalar host payloads must already be registered with the evaluator
template <typename T>
ValueTypeInfo hostResultType(const Evaluator& evaluator)
{
  using Result=
    std::remove_cvref_t<T>;
  ValueTypeInfo result;
  if constexpr(std::is_void_v<Result>)
    result.cppType=typeid(std::monostate);
  else if constexpr(std::is_same_v<Result,bool>)
    result.cppType=typeid(int);
  else if constexpr(std::is_same_v<Result,int> or std::is_same_v<Result,int64_t> or
                    std::is_same_v<Result,double> or std::is_same_v<Result,std::string>)
    result.cppType=typeid(Result);
  else if constexpr(not std::is_constructible_v<ValueT,Result> and not std::is_same_v<Result,Value>)
    {
      if(const auto found=evaluator.hostObjectTypes.find(typeid(Result));found!=evaluator.hostObjectTypes.end())
        result={typeid(Result),found->second};
    }
  
  return result;
}

/// Install a typed C++ lambda as a script function
/// The provided names correspond to its parameters and enable named arguments
/// An optional leading Evaluator& receives the caller context without a script argument
/// An optional description is displayed by verbose entity discovery
/// Native return types are retained automatically; erased Value/HObject returns may supply resultType
/// Its captured references must outlive the evaluator's registered function
template <typename Callable>
void registerHostFunction(Evaluator& evaluator,
                          const std::string& name,
                          std::vector<std::string> argumentNames,
                          Callable callable,
                          std::string description={},
                          ValueTypeInfo resultType={})
{
  using Signature=
    HostCallableSignature<decltype(&Callable::operator())>;
  
  if(name.empty())
    errorEmitter("host function name must not be empty");
  
  if(argumentNames.size()!=Signature::arity)
    errorEmitter("argument name count does not match signature of ",name);
  
  for(size_t i=0;i<argumentNames.size();i++)
    {
      if(argumentNames[i].empty())
        errorEmitter("empty argument name in ",name);
      
      for(size_t j=0;j<i;j++)
        if(argumentNames[i]==argumentNames[j])
          errorEmitter("duplicate argument name ",argumentNames[i]," in ",name);
    }
  
  if(not resultType.cppType)
    resultType=hostResultType<typename Signature::Result>(evaluator);
  
  (*evaluator.curEnv)[name]=HostFunction{
    [callable=std::move(callable),name](std::vector<std::shared_ptr<Value>>& values,
                                      Evaluator& caller) mutable -> Value
    {
      return Signature::invoke(callable,values,name.c_str(),caller);
    },std::move(argumentNames),false,std::move(description),
    Signature::arity==0?FunctionArgumentForms(1):FunctionArgumentForms{},std::move(resultType)};
}

/// Registers Type under its C++ identifier, using its static members() function
/// Use an unqualified type name or alias that is also a valid script identifier
/// Returns the shared metadata, just like registerHObject
#define REGISTER_HOBJECT(EV,TYPE) \
  ::registerHObject<TYPE>((EV),#TYPE,TYPE::members())

#endif
