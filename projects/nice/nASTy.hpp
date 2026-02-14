#include <cmath>

#include <parsePact.hpp>

using namespace pp::internal;

inline auto getParser()
{
#define BARE_WHITESPACES			\
  "( |\\n|\\t)+"
  
#define CPP_COMMENT				\
  "//[^\\n]+\\n"
  
#define C_COMMENT				\
  "/\\*([^\\*]|\\*[^/])+\\*/"
    
    const char cGrammar[]=
    "c {"
    "   %whitespace \"" BARE_WHITESPACES "|" CPP_COMMENT "|" C_COMMENT "\";"
    "   %none lowerThanElse;"
    "   %none \"else\";"
    "   %left \"\\(\";"
    "   %left \"\\)\";"
    "   %none BRACKETS;"
    "   %left \",\";"
    "   %right \"=\";"
    "   %right \"\\+=\";"
    "   %right \"\\-=\";"
    "   %right \"\\*=\";"
    "   %right \"/=\";"
    "   %left \"\\|\\||or\";"
    "   %left \"&&|and\";"
    "   %left \"==\";"
    "   %left \"!=\";"
    "   %left \"<\";"
    "   %left \"<=\";"
    "   %left \">\";"
    "   %left \">=\";"
    "   %left \"\\+\";"
    "   %left \"\\-\";"
    "   %left \"%\";"
    "   %left \"\\*\";"
    "   %left \"/\";"
    "   %right UNARY_ARITHMETIC;"
    "   %right \"!\";"
    "   %left \"\\-\\-\";"
    "   %left \"\\+\\+\";"
    "   %none FUNCTION_CALL;"
    ""
    "   statement : expression_statement [return($0)]"
    "             | compound_statement [return($0)]"
    "             | forStatement [return($0)]"
    "             | ifStatement [return($0)]"
    "             | functionDefinition [return($0)]"
    "             | structDefinition [return($0)]"
    "             | \"return\" expression_statement [funcReturn($1)]"
    "             ;"
    "   functionDefinition : \"fun\" identifier \"\\(\" variadicFuncParList \"\\)\" compound_statement [funcDef($1,$3,$5)]"
    "                      ;"
    "   structDefinition : \"struct\" identifier \"{\" structMemberList \"}\" \";\" [structDef($1,$3)]"
    "                    ;"
    "   structMemberList : [emptyStructMemberList]"
    "                    | nonEmptyStructMemberList [return($0)]"
    "                    ;"
    "   nonEmptyStructMemberList : structMember \";\" [firstStructMemberOfList($0)]"
    "                            | nonEmptyStructMemberList structMember \";\" [appendStructMemberToList($0,$1)]"
    "                            ;"
    "   structMember : identifier [structMemberCreate($0)]"
    "                | identifier \"=\" expression [structMemberCreateWithDef($0,$2)]"
    "                ;"
    "   variadicFuncParList: funcParList [return($0)]"
    "                      | variadicFuncPar [emptyVariadicParList($0)]"
    "                      | nonEmptyFuncParList \",\" variadicFuncPar [makeVariadicParList($0)]"
    "                      ;"
    "   variadicFuncPar: \"\\.\\.\\.\" [return(1)]"
    "                  | \"&\" \"\\.\\.\\.\" [return(2)]"
    "                  ;"
    "   funcParList: [emptyFuncParList]"
    "              | nonEmptyFuncParList [return($0)]"
    "              ;"
    "   nonEmptyFuncParList: funcPar [firstFuncParOfList]"
    "                      | nonEmptyFuncParList \",\" funcPar [appendFuncParToList($0,$2)]"
    "                       ;"
    "   funcPar: neededFuncPar [return($0)]"
    "          | neededFuncPar \"=\" expression [addFuncParDefault($0,$2)]"
    "          ;"
    "   neededFuncPar: identifier [funcParCreate(0,$0)]"
    "                | \"&\" identifier [funcParCreate(1,$1)]"
    "                ;"
    "   forStatement: \"for\" \"\\(\" forInit \";\" forCheck \";\" forIncr \"\\)\" statement [forStatement($2,$4,$6,$8)]"
    "               ;"
    "   forInit: expression [return($0)]"
    "          |"
    "          ;"
    "   forCheck: expression [return($0)]"
    "           |"
    "           ;"
    "   forIncr: expression [return($0)]"
    "          |"
    "          ;"
    "   ifStatement: \"if\" \"\\(\" expression \"\\)\" statement %precedence lowerThanElse [ifStatement($2,$4)]"
    "              | \"if\" \"\\(\" expression \"\\)\" statement \"else\" statement [ifElseStatement($2,$4,$6)]"
    "              ;"
    "    compound_statement : \"{\" statements \"}\" [return($1)]"
    "                       ;"
    "    statements : [createStatements]"
    "               | statements statement [appendStatement]"
    "               ;"
    "    expression_statement : expression \";\" [return($0)]"
    "                         ;"
    "    expression : postfix_expression %precedence \"\\+\\+\" [return($0)]"
    "               | expression \"\\+\" expression [binarySum($0,$2)]"
    "               | expression \"\\-\" expression [binaryDiff($0,$2)]"
    "               | expression \"\\*\" expression [binaryProd($0,$2)]"
    "               | expression \"/\" expression [binaryDiv($0,$2)]"
    "               | expression \"%\" expression [binaryMod($0,$2)]"
    "               | \"\\+\" expression %precedence UNARY_ARITHMETIC [unaryPlus($1)]"
    "               | \"\\-\" expression %precedence UNARY_ARITHMETIC [unaryMinus($1)]"
    "               | \"!\" expression [unaryNot($1)]"
    "               | \"&\" expression %precedence UNARY_ARITHMETIC [unaryReference($1)]"
    "               | expression \"<\" expression [binarySmaller($0,$2)]"
    "               | expression \">\" expression [binaryGreater($0,$2)]"
    "               | expression \"<=\" expression [binarySmallerEqual($0,$2)]"
    "               | expression \">=\" expression [binaryGreaterEqual($0,$2)]"
    "               | expression \"==\" expression [binaryCompare($0,$2)]"
    "               | expression \"!=\" expression [binaryInequal($0,$2)]"
    "               | expression \"&&|and\" expression [binaryAnd($0,$2)]"
    "               | expression \"\\|\\||or\" expression [binaryOr($0,$2)]"
    "               | assign_expression [return($0)]"
    "               ;"
    "    postfix_expression : primary_expression [return($0)]"
    "                       | postfix_expression \"\\+\\+\" [unaryPostfixIncrement($0)]"
    "                       | postfix_expression \"\\-\\-\" [unaryPostfixDecrement($0)]"
    "                       | postfix_expression \"\\(\" funcArgsList \"\\)\" %precedence FUNCTION_CALL [funcCall($0,$2)] "
    "                       | postfix_expression \"\\[\" expression \"\\]\" [subscribe($0,$2)] "
    "                       | postfix_expression \"\\.\" identifier [memberAccess($0,$2)] "
    "                       ;"
    "    funcArgsList : [emptyFuncArgList]"
    "                 | nonemptyFuncArgsList [return($0)]"
    "                 ;"
    "    nonemptyFuncArgsList : arg [firstFuncArgOfList]"
    "                         | funcArgsList \",\" arg [appendFuncArgToList($0,$2)]"
    "                         ;"
    "    arg : \"\\.\" identifier \"=\" expression [createFuncArg($1,$3)]"
    "        | expression [createFuncArg(\"\",$0)]"
    "        ;"
    "    primary_expression : identifier [return($0)]"
    "                       | \"[0-9]+\" [convToInt($0)]"
    "                       | \"([0-9]+(\\.[0-9]*)?|(\\.[0-9]+))((e|E)(\\+|\\-)?[0-9]+)?\" [convToFloat($0)]"
    "                       | \"\\\"[^\\\"]*\\\"\" [convToStr($0)]"
    "                       | \"\\(\" expression \"\\)\" %precedence BRACKETS [return($1)]"
    "                       | \"lambda\" \"\\(\" variadicFuncParList \"\\)\" compound_statement [lambdaFuncDef($2,$4)]"
    "                       ;"
    "    assign_expression : postfix_expression \"=\" expression %precedence \"=\" [assign($0,$2,0)]"
    "                      | postfix_expression \"\\*=\" expression [assign($0,$2,1)]"
    "                      | postfix_expression \"/=\" expression [assign($0,$2,2)]"
    "                      | postfix_expression \"\\+=\" expression [assign($0,$2,3)]"
    "                      | postfix_expression \"\\-=\" expression [assign($0,$2,4)]"
    "                      ;"
    "    identifier : \"[a-zA-Z_][a-zA-Z0-9_]*\" [convToId($0)]"
    "               ;"
    "}";
  
  return
    pp::createGrammar(cGrammar);
}

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

struct ValuesList;

struct StructDef;

struct Struct;

struct FuncArg;

using Value=
  std::variant<std::monostate,std::string,int,double,Function,FuncArg,HostFunction,ValueRef,ValuesList,StructDef,Struct>;

struct HostedRef
{
  std::function<Value()> get;
  
  std::function<void(const AssignMode&,const Value&)> set;
};

struct ValueRef
{
  std::shared_ptr<Value> ref;
};

struct ValuesList
{
  std::vector<std::shared_ptr<Value>> data;
  
  Value& operator[](const int& i);
};

struct Function
{
  const FuncNode* funcNode;
  
  std::shared_ptr<Environment> env;
};

struct HostFunction :
  std::function<Value(std::vector<std::shared_ptr<Value>>&,Evaluator&)>
{
};

struct StructDef
{
  const StructDefNode* structDefNode;
};

struct Struct
{
  const StructDefNode* structDefNode;
  
  std::vector<std::shared_ptr<Value>> members;
};

struct IdNode
{
  std::string name;
};

struct FuncArgNode
{
  std::string name;
  
  std::shared_ptr<ASTNode> expr;
};

struct FuncArg
{
  std::string name;
  
  std::shared_ptr<Value> value;
};

struct FuncArgListNode
{
  std::vector<FuncArgNode> list;
};

struct SubscribeNode
{
  std::shared_ptr<ASTNode> base;
  
  std::shared_ptr<ASTNode> subscr;
};

struct MemberAccessNode
{
  std::shared_ptr<ASTNode> base;
  
  std::string member;
};

struct ValueNode
{
  Value value;
};

struct ASTNodesNode
{
  std::vector<std::shared_ptr<ASTNode>> list;
};

struct IfNode
{
  std::vector<std::shared_ptr<ASTNode>> subNodes;
};

struct FuncParNode
{
  std::string name;
  
  bool isRef;
  
  std::shared_ptr<ASTNode> def;
  
  static constexpr const char VARIADIC_PAR_NAME[]="__VA__ARGS__";
  
  static FuncParNode getVariadic(const bool& isRef)
  {
    return {.name=VARIADIC_PAR_NAME,.isRef=isRef,.def=std::make_shared<ASTNode>(ValueNode{ValuesList{}})};
  }
};

template <typename D>
struct Searchable
{
  auto iParOfName(const std::string& name) const
  {
    const auto& list=((const D*)this)->list;
    
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

struct FuncParListNode:
  Searchable<FuncParListNode>
{
  std::vector<FuncParNode> list;
  
  bool isVariadic() const
  {
    return list.size() and list.back().name==FuncParNode::VARIADIC_PAR_NAME;
  }
  
  void makeVariadic(const int& isRef)
  {
    if(isVariadic())
      errorEmitter("Cannot make a list of pars variadic twice");
    
    list.push_back(FuncParNode::getVariadic(isRef));
    
  };
  
  static FuncParListNode getEmptyVariadic(const int& isRef)
  {
    FuncParListNode list;
    
    list.makeVariadic(isRef);
    
    return list;
  }
};

struct ForNode
{
  std::vector<std::shared_ptr<ASTNode>> subNodes;
};

struct FuncNode
{
  FuncParListNode pars;
  
  std::shared_ptr<ASTNode> body;
};

struct FuncDefNode
{
  std::string name;
  
  std::shared_ptr<FuncNode> fun;
};

struct FuncCallNode
{
  std::string fun;
  
  FuncArgListNode args;
};

struct StructMemberNode
{
  std::string name;
  
  std::shared_ptr<ASTNode> def;
};

struct StructMemberListNode:
  Searchable<StructMemberListNode>
{
  std::vector<StructMemberNode> list;
};

struct StructDefNode
{
  std::string name;
  
  StructMemberListNode members;
};

struct AssignNode
{
  std::shared_ptr<ASTNode> lhs;
  
  std::shared_ptr<ASTNode> rhs;
  
  AssignMode mode;
};

struct PrePostfixOpNode
{
  std::shared_ptr<ASTNode> arg;
  
  std::function<Value(Value&)> op;
};

struct UnaryRefNode
{
  std::shared_ptr<ASTNode> arg;
};

struct UnOpNode
{
  std::shared_ptr<ASTNode> arg;
  
  std::function<Value(const Value&)> op;
};

struct BinOpNode
{
  std::shared_ptr<ASTNode> arg1;
  
  std::shared_ptr<ASTNode> arg2;
  
  std::function<Value(const Value&,const Value&)> op;
};

struct ReturnNode
{
  std::shared_ptr<ASTNode> arg;
};

struct Environment
{
  std::shared_ptr<Environment> parent;
  
  std::unordered_map<std::string,std::shared_ptr<Value>> varTable;
  
  Environment& uppermost()
  {
    Environment* res=this;
    
    while(Environment* par=&*res->parent)
      res=par;
    
    return *res;
  }
  
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
  
  Value& at(const std::string& name)
  {
    std::shared_ptr<Value> f=find(name);
    
    if(f==nullptr)
      errorEmitter("Trying to use uninitialized variable ",name);
    
    return *f;
  }
  
  std::pair<std::shared_ptr<Value>,bool> getRefOrInsert(const std::string& name)
  {
    if(auto f=find(name))
      return {f,false};
    else
      return {varTable[name]=std::make_shared<Value>(std::monostate{}),true};
  }
  
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

struct CallFrame
{
  std::shared_ptr<Environment> env;
  
  CallFrame(std::shared_ptr<Environment> parentEnv=nullptr) :
    env{std::make_shared<Environment>(parentEnv)}
  {
  }
};

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

template <typename T>
T& fetch(std::vector<std::shared_ptr<ASTNode>>& list,
	 const size_t& i)
{
  using namespace pp::internal;
  
  if(const size_t n=list.size();n<i)
    errorEmitter(n," nodes received, aksed for node #",i);
  
  return fetch<T>(list[i],std::to_string(i).c_str());
}

template <typename T,
	  typename L>
std::shared_ptr<ASTNode> makeFirstOfList(std::shared_ptr<ASTNode>& elem)
{
  return std::make_shared<ASTNode>(L{.list{fetch<T>(elem)}});
}

template <typename T,
	  typename L>
std::shared_ptr<ASTNode>& appendToList(std::shared_ptr<ASTNode>& list,
				       std::shared_ptr<ASTNode>& elem)
{
 fetch<L>(list).list.push_back(std::move(fetch<T>(elem)));
 
 return list;
}

template <typename L>
std::shared_ptr<ASTNode>& mergeLists(std::shared_ptr<ASTNode>& list1,
				     std::shared_ptr<ASTNode>& list2)
{
  auto& first=fetch<L>(list1).list;
  auto& second=fetch<L>(list2).list;
  
  first.insert(first.end(),std::make_move_iterator(second.begin()),std::make_move_iterator(second.end()));
  
  return list1;
}

inline auto getParseTreeExecuctor(const std::vector<std::string_view>& requiredActions)
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
    [](std::vector<std::shared_ptr<ASTNode>>& subNodes) -> std::shared_ptr<ASTNode> \
    {									\
      ENSURE_N_SYMBOLS(NAME,N);						\
      									\
      BODY;								\
    }
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("createStatements",0,return std::make_shared<ASTNode>(ASTNodesNode{}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("appendStatement",2,
				if(ASTNodesNode* l=std::get_if<ASTNodesNode>(&*(subNodes[0]));l==nullptr)
				  errorEmitter("first argument is not a list of statement: ",
					       std::visit([](const auto& x)
					       {
						 return typeid(decltype(x)).name();
					       },*(subNodes[0])));
				else
				  l->list.emplace_back(subNodes[1]);
				return subNodes[0]);
  PROVIDE_ACTION_WITH_N_SYMBOLS("convToInt",1,return std::make_shared<ASTNode>(ValueNode{atoi(unvariant<std::string>(fetch<ValueNode>(subNodes,0).value).c_str())}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("convToStr",1,return std::make_shared<ASTNode>(ValueNode{unescapeString(unvariant<std::string>(fetch<ValueNode>(subNodes,0).value))}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("convToId",1,return std::make_shared<ASTNode>(IdNode{.name=unvariant<std::string>(fetch<ValueNode>(subNodes,0).value)}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("convToFloat",1,return std::make_shared<ASTNode>(ValueNode{strtod(unvariant<std::string>(fetch<ValueNode>(subNodes,0).value).c_str(),nullptr)}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("return",1,return subNodes[0]);
  PROVIDE_ACTION_WITH_N_SYMBOLS("createFuncArg",2,
				std::string name{};
				if(IdNode* id=std::get_if<IdNode>(&*subNodes[0]))
				  name=id->name;
				return std::make_shared<ASTNode>(FuncArgNode{.name=name,.expr=subNodes[1]}));
  
#define DEFINE_UNOP(OP,NAME)						\
  PROVIDE_ACTION_WITH_N_SYMBOLS("unary" #NAME,1,return std::make_shared<ASTNode>(UnOpNode{.arg=subNodes[0], \
	  .op=[](const Value& arg)					\
	  {								\
	    return std::visit([](const auto& arg) -> Value		\
	    {								\
	      if constexpr(requires {OP arg;})				\
		return OP arg;						\
	      else							\
		{							\
		  errorEmitter("Cannot " #NAME" the type: %s",typeid(arg).name()); \
		  return std::monostate{};				\
		}},							\
			      arg);				\
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
		return arg1 OP arg2;					\
	      else							\
		{							\
		  errorEmitter("Cannot " #NAME" the types: %s %s",typeid(arg1).name(),typeid(arg2).name()); \
		  return std::monostate{};				\
		}},							\
			      arg1,arg2);				\
	  }}))
  
  DEFINE_BINOP(+,Sum);
  DEFINE_BINOP(-,Diff);
  DEFINE_BINOP(-,Prod);
  DEFINE_BINOP(/,Div);
  DEFINE_BINOP(%,Mod);
  DEFINE_BINOP(<,Smaller);
  DEFINE_BINOP(>,Greater);
  DEFINE_BINOP(<=,SmallerEqual);
  DEFINE_BINOP(>=,GreaterEqual);
  DEFINE_BINOP(==,Compare);
  DEFINE_BINOP(!=,Inequal);
  DEFINE_BINOP(or,Or);
  DEFINE_BINOP(and,And);
  
#undef DEFINE_BINOP
  
#define DEFINE_PREPOSTFIX_OP(NAME,PRE,POST)					\
  PROVIDE_ACTION_WITH_N_SYMBOLS("unary" #NAME,1,return std::make_shared<ASTNode>(PrePostfixOpNode{.arg=subNodes[0], \
	  .op=[](Value& arg)						\
	  {								\
	    return std::visit([](auto& arg) -> Value			\
	    {								\
	      if constexpr(requires {PRE arg POST;})			\
		return PRE arg POST;					\
	      else							\
		{							\
		  errorEmitter("Cannot " #NAME" the type: %s",typeid(arg).name()); \
		  return std::monostate{};				\
		}},							\
			      arg);				\
	  }}))
  
  DEFINE_PREPOSTFIX_OP(PostfixIncrement,,++);
  //  DEFINE_PREPOSTFIX_OP(PrefixIncrement,++,);
  DEFINE_PREPOSTFIX_OP(PostfixDecrement,,--);
  //DEFINE_PREPOSTFIX_OP(PrefixDecrement,--,);
  
#undef DEFINE_PREPOSTFIX_OP
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("structMemberCreate",1,return std::make_shared<ASTNode>(StructMemberNode{.name=fetch<IdNode>(subNodes,0).name}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("structMemberCreateWithDef",2,return std::make_shared<ASTNode>(StructMemberNode{.name=fetch<IdNode>(subNodes,0).name,
														.def=subNodes[1]}));
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("funcParCreate",2,return std::make_shared<ASTNode>(FuncParNode{.name=fetch<IdNode>(subNodes,1).name,
											   .isRef=(bool)unvariant<int>(fetch<ValueNode>(subNodes,0).value)}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("addFuncParDefault",2,const FuncParNode& in=fetch<FuncParNode>(subNodes,0);
				return std::make_shared<ASTNode>(FuncParNode{.name=in.name,.isRef=in.isRef,.def=subNodes[1]}));
#define PROVIDE_FUNC_LIST_ACTIONS(NAME)					\
  PROVIDE_ACTION_WITH_N_SYMBOLS("empty" #NAME "List",0,return std::make_shared<ASTNode>(NAME ## ListNode{})); \
  PROVIDE_ACTION_WITH_N_SYMBOLS("first" #NAME "OfList",1,return makeFirstOfList<NAME ## Node,NAME ## ListNode>(subNodes[0])); \
  PROVIDE_ACTION_WITH_N_SYMBOLS("append" #NAME "ToList",2,return appendToList<NAME ## Node,NAME ## ListNode>(subNodes[0],subNodes[1])); \
  
  PROVIDE_FUNC_LIST_ACTIONS(FuncPar);
  PROVIDE_FUNC_LIST_ACTIONS(FuncArg);
  PROVIDE_FUNC_LIST_ACTIONS(StructMember);
  
#undef PROVIDE_FUNC_LIST_ACTIONS
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("emptyVariadicParList",1,return std::make_shared<ASTNode>(FuncParListNode::getEmptyVariadic(unvariant<int>(fetch<ValueNode>(subNodes,0).value))));
  PROVIDE_ACTION_WITH_N_SYMBOLS("makeVariadicParList",2,
				fetch<FuncParListNode>(subNodes,0).makeVariadic(unvariant<int>(fetch<ValueNode>(subNodes,1).value));
				return subNodes[0]);
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("funcDef",3,
				return std::make_shared<ASTNode>(FuncDefNode{.name=fetch<IdNode>(subNodes,0).name,
									     .fun=std::make_shared<FuncNode>(fetch<FuncParListNode>(subNodes,1),
													     subNodes[2])}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("structDef",2,
				return std::make_shared<ASTNode>(StructDefNode{.name=fetch<IdNode>(subNodes,0).name,
									     .members=fetch<StructMemberListNode>(subNodes,1)}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("lambdaFuncDef",2,
				return std::make_shared<ASTNode>(FuncNode{fetch<FuncParListNode>(subNodes,0),
									  subNodes[1]}));
  
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("funcCall",2,return std::make_shared<ASTNode>(FuncCallNode{.fun=fetch<IdNode>(subNodes,0).name,
											   .args=fetch<FuncArgListNode>(subNodes,1)}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("funcReturn",1,return std::make_shared<ASTNode>(ReturnNode{.arg=subNodes[0]}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("subscribe",2,return std::make_shared<ASTNode>(SubscribeNode{.base=subNodes[0],.subscr=subNodes[1]}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("memberAccess",2,return std::make_shared<ASTNode>(MemberAccessNode{.base=subNodes[0],.member=fetch<IdNode>(subNodes,1).name}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("forStatement",4,return std::make_shared<ASTNode>(ForNode{.subNodes{subNodes}}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("ifStatement",2,return std::make_shared<ASTNode>(IfNode{.subNodes{subNodes}}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("ifElseStatement",3,return std::make_shared<ASTNode>(IfNode{.subNodes{subNodes}}));
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("assign",3,return std::make_shared<ASTNode>(AssignNode{.lhs=subNodes[0],
	  .rhs=subNodes[1],
	  .mode=(AssignMode)unvariant<int>(fetch<ValueNode>(subNodes,2).value)}));
  
#undef PROVIDE_ACTION_WITH_N_SYMBOLS
  
  return ParseTreeExecutor(providedActions,requiredActions);
}

template <typename L,
	  typename R>
void trySet(L& lhs,
	    const AssignMode& mode,
	    const R& rhs)
{
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

struct Evaluator
{
  struct EvalResult
  {
    using Ref=
      std::variant<std::monostate,ValueRef,HostedRef>;
    
    Ref ref;
    
    bool hasRef() const
    {
      return not std::holds_alternative<std::monostate>(ref);
    }
    
    std::optional<Value> value;
    
    void set(const AssignMode& mode,
	     const Value& v)
    {
      std::visit(Overload{[&v,
			   &mode](const ValueRef& vr)
      {
	trySet(*vr.ref,mode,v);
      },
	    [&v,
	     &mode](const HostedRef& hr)
	    {
	      hr.set(mode,v);
	    },
	    [](const std::monostate&)
	    {
	      errorEmitter("trying to set a non lhs");
	    }},ref);
    }
    
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
    
    static void ensureInitialized(const Value& v)
    {
      if(std::get_if<std::monostate>(&v))
	errorEmitter("using uninitialized value");
    }
    
    Value asRhs() const
    {
      if(hasRef())
	return
	  std::visit(Overload{[this](const ValueRef& vr)->Value
	  {
	    if(const std::shared_ptr<Value>& r=vr.ref)
	      return *r;
	    else
	    {
	      if(not value)
		errorEmitter("trying to acces an empty evaluator result");
	      
	      ensureInitialized(*value);
	      
	      return *value;
	    }
	  },
		[](const HostedRef& hr)->Value
		{
		  return hr.get();
		},
		[](const std::monostate&) ->Value
		{
		  errorEmitter("impossible");
		  
		  return {};
		}},ref);
      else
	return *value;
    }
  };
  
  std::vector<CallFrame> callFrames;
  
  std::shared_ptr<Environment> curEnv;
  
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
  
  EnvironmentWindUnwinder pushNewEnvironment()
  {
    return {*this,std::make_shared<Environment>(curEnv)};
  }
  
  Evaluator()
  {
    curEnv=callFrames.emplace_back().env;
    
    Environment& env=*curEnv;
    
    env["M_PI"]=M_PI;
    
    env["print"]=
      HostFunction{
      [](std::vector<std::shared_ptr<Value>>& args,
	 Evaluator&)->Value
      {
	for(std::shared_ptr<Value>& arg : args)
	  std::visit([](const auto& v)
	  {
	    if constexpr(requires {std::cout<<v;})
	      std::cout<<v;
	    else
	      std::cout<<"unprintable type: "<<typeid(decltype(v)).name()<<"\n";
	  },*arg);
	
	return {};
      }};
    
    env["makeVec"]=
      HostFunction{
      [](std::vector<std::shared_ptr<Value>>& args,
	 Evaluator&)->Value
      {
	return ValuesList{args};
      }};
    
    env["forEachEl"]=
      HostFunction{
      [](std::vector<std::shared_ptr<Value>>& args,
	 Evaluator& ev)->Value
      {
	if(const size_t n=args.size();n!=2)
	  errorEmitter("trying to call forEachEl with ",n," elements in place of 2");
	
	std::vector<std::shared_ptr<Value>> v(1);
	if(ValuesList* list=std::get_if<ValuesList>(&*args[0]))
	  for(std::shared_ptr<Value> val : list->data)
	    {
	      v[0]=val;
	      
	      ev.callFunction(*args[1],v);
	    }
	else
	  errorEmitter("trying to call \"forEachEl\" on a non-vector type ",variantInnerTypeName(*args[0]));
	
	return {};
      }};
    
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
		int* r=std::get_if<int>(&*args[i]);
		if(not r)
		  errorEmitter("arg ",i," is of type ",variantInnerTypeName(*args[i])," instead of ",typeid(int).name()," (int)");
		
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
	
	return v;
      }};
    
#define REGISTER_ARGLESS_HOST_FUNCTION(NAME)		\
    env[#NAME]=						\
      HostFunction{					\
      [](std::vector<std::shared_ptr<Value>>& args,\
	 Evaluator&)->Value			   \
      {						   \
	return NAME();				   \
      }}
    
#define REGISTER_HOST_FUNCTION(NAME,ARGS...)				\
    env[#NAME]=								\
      HostFunction{							\
      [](std::vector<std::shared_ptr<Value>>& args,			\
	 Evaluator&)->Value						\
      {									\
	constexpr size_t N=						\
	  N_VARIADIC_ARGS(ARGS);					\
    									\
	const size_t n=							\
	  args.size();							\
									\
	if(N!=n)							\
	  errorEmitter("trying to call function ",#NAME,		\
		       " which expects ",N," args with ",n);		\
									\
	return std::visit([](auto&...args) ->Value			\
	{								\
	  if constexpr(requires {NAME(std::forward<decltype(args)>(args)...);}) \
	    if constexpr(std::is_same_v<void,decltype(NAME(std::forward<decltype(args)>(args)...))>) \
	      NAME(std::forward<decltype(args)>(args)...);		\
	    else							\
	      return NAME(std::forward<decltype(args)>(args)...);	\
	  else								\
	    errorEmitter("Trying to call ",#NAME,		\
				       " function with impossible args ",typeid(args).name()...); \
	  								\
	  return std::monostate{};					\
	}								\
	  ,ARGS);							\
      }}
    
    REGISTER_ARGLESS_HOST_FUNCTION(rand);
    
    REGISTER_HOST_FUNCTION(srand,*args[0]);
    REGISTER_HOST_FUNCTION(sqrt,*args[0]);
    REGISTER_HOST_FUNCTION(exp,*args[0]);
    REGISTER_HOST_FUNCTION(sin,*args[0]);
    REGISTER_HOST_FUNCTION(cos,*args[0]);
    REGISTER_HOST_FUNCTION(tan,*args[0]);
    REGISTER_HOST_FUNCTION(pow,*args[0],*args[1]);
    
    {
      using namespace std;
      REGISTER_HOST_FUNCTION(to_string,*args[0]);
    }
  }
  
  EvalResult createStruct(const StructDef& def,
			  std::vector<std::shared_ptr<Value>>& args)
  {
    const std::string& name=def.structDefNode->name;
    
    diagnostic("Creating struct ",name,"\n");
    
    const size_t nArgs=args.size();
    const size_t nMembers=def.structDefNode->members.list.size();
    
    if(args.size()!=nMembers)
      errorEmitter("Calling struct ",name," constructor with ",nArgs," args when expecting ",nMembers);
      
    Struct res{.structDefNode=def.structDefNode,.members=args};
    
    return {.value=res};
  }
  
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
	    if(ValueRef* vr=std::get_if<ValueRef>(&*arg))
	      arg=vr->ref;
	    else
	      errorEmitter("parameter ",par.name," has to be taken by ref but passed argument is not a ref");
	  }
	
	if(not curEnv->varTable.try_emplace(par.name,arg).second)
	  errorEmitter("argument ",par.name," already specified");
      }
    
    return std::visit(*this,*fn.body);
  }
  
  EvalResult callFunction(const Value& vFun,
			  std::vector<std::shared_ptr<Value>>& args)
  {
    if(const Function* f=std::get_if<Function>(&vFun))
      return callEmbeddedFunction(*f,args);
    else
      if(const HostFunction* hf=std::get_if<HostFunction>(&vFun))
	return {.value=(*hf)(args,*this)};
      else
	if(const StructDef* sd=std::get_if<StructDef>(&vFun))
	  return createStruct(*sd,args);
	else
	  errorEmitter("Variable is of type ",variantInnerTypeName(vFun)," not supported by the callFunction backend");
    
    return {};
  }
  
  EvalResult callFunction(const std::string& name,
			  std::vector<std::shared_ptr<Value>>& args)
  {
    diagnostic("Going to call function ",name.c_str(),"\n");
    
    return callFunction(curEnv->at(name),args);
  }
  
  EvalResult operator()(const ValueNode& valueNode)
  {
    return {.value=valueNode.value};
  }
  
  EvalResult operator()(const IdNode& idNode)
  {
    return {.ref=ValueRef(curEnv->getRefOrInsert(idNode.name).first)};
  }
  
  EvalResult operator()(const UnaryRefNode& uRNode)
  {
    EvalResult tmp=std::visit(*this,*uRNode.arg);
    
    ValueRef* _vr=std::get_if<ValueRef>(&tmp.ref);
    if(not _vr)
      errorEmitter("cannot take the reference of a non-embedded object");
    
    ValueRef& vr=*_vr;
    if(vr.ref)
      errorEmitter("subNode has not returned a reference");
      
    return {.value{vr}};
  }
  
  EvalResult operator()(const StructDefNode& structDefNode)
  {
    const std::string& name=
      structDefNode.name;
    
    if(curEnv->find(name))
      errorEmitter("Redefining a struct which has a name ",name," already defined");
    
    (*curEnv)[name]=StructDef{.structDefNode=&structDefNode};
    
    return {};
  }
  
  EvalResult operator()(const SubscribeNode& subscribeNode)
  {
    EvalResult _b=std::visit(*this,*subscribeNode.base);
    const Value s=std::visit(*this,*subscribeNode.subscr).asRhs();
    Value& b=_b.hasRef()?_b.asLhs():*_b.value;
    
    if(ValuesList* vl=std::get_if<ValuesList>(&b))
      if(const int* mi=std::get_if<int>(&s))
	{
	  int i=*mi;
	  
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
	errorEmitter("index of subscrition is not an integer but is of type: ",variantInnerTypeName(s));
    else
      errorEmitter("trying to subscribe non-list of type: ",variantInnerTypeName(b));
    
    return {};
  }
  
  EvalResult operator()(const MemberAccessNode& memberAccessNode)
  {
    EvalResult _b=std::visit(*this,*memberAccessNode.base);
    Value& b=_b.hasRef()?_b.asLhs():*_b.value;
    
    const std::string& memberName=memberAccessNode.member;
    
    if(Struct* es=std::get_if<Struct>(&b))
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
      errorEmitter("trying to access member ",memberName," of a non-stuct of type: ",variantInnerTypeName(b));
    
    return {};
  }
  
  EvalResult operator()(const ForNode& forNode)
  {
    EnvironmentWindUnwinder windUnwind=pushNewEnvironment();
    
    for(std::visit(*this,*forNode.subNodes[0]);
    	std::visit([](const auto& v)
	{
	  if constexpr(std::is_convertible_v<decltype(v),bool>)
	    return (bool)v;
	  else
	    errorEmitter("Cannot convert the type to bool");
	  
	  return false;
	},std::visit(*this,*forNode.subNodes[1]).asRhs());
	std::visit(*this,*forNode.subNodes[2]).asRhs())
      std::visit(*this,*forNode.subNodes[3]);
    
    return {};
  }
  
  EvalResult operator()(const IfNode& ifNode)
  {
    EnvironmentWindUnwinder windUnwind=pushNewEnvironment();
    
    if(std::visit([](const auto& v)
	{
	  if constexpr(std::is_convertible_v<decltype(v),bool>)
	    return (bool)v;
	  else
	    errorEmitter("Cannot convert the type to bool");
	  
	  return false;
	},std::visit(*this,*ifNode.subNodes[0]).asRhs()))
      std::visit(*this,*ifNode.subNodes[1]);
    else
      if(ifNode.subNodes.size()>=2)
	std::visit(*this,*ifNode.subNodes[2]);
    
    return {};
  }
  
  EvalResult operator()(const FuncParNode&)
  {
    return {};
  }
  
  EvalResult operator()(const FuncParListNode&)
  {
    return {};
  }
  
  EvalResult operator()(const StructMemberNode& structMemberNode)
  {
    if(std::shared_ptr<ASTNode> def=structMemberNode.def)
      return {.value=std::visit(*this,*def).asRhs()};
    else
      return {.value=std::monostate{}};
  }
  
  EvalResult operator()(const StructMemberListNode& structMemberListNode)
  {
    ValuesList memberList;
    
    const size_t nMembers=structMemberListNode.list.size();
    
    memberList.data.reserve(nMembers);
    
    diagnostic("Preparing a list of ",nMembers," members\n");
    
    for(const StructMemberNode& structMemberNode : structMemberListNode.list)
      memberList.data.push_back(std::make_shared<Value>((*this)(structMemberNode).asRhs()));
    
    return {.value=memberList};
  }
  
  EvalResult operator()(const FuncArgNode& argNode)
  {
    return {.value=FuncArg{.name=argNode.name,
			   .value=std::make_shared<Value>(std::visit(*this,*argNode.expr).asRhs())}};
  }
  
  EvalResult operator()(const FuncArgListNode& argListNode)
  {
    ValuesList argList;
    
    const size_t nArgs=argListNode.list.size();
    
    argList.data.reserve(nArgs);
    
    diagnostic("Preparing a list of ",nArgs," args\n");
    
    for(const FuncArgNode& argNode : argListNode.list)
      argList.data.push_back(std::make_shared<Value>((*this)(argNode).asRhs()));
    
    return {.value=argList};
  }
  
  EvalResult operator()(const FuncCallNode& funcCallNode)
  {
    Value __args=(*this)(funcCallNode.args).asRhs();
    
    ValuesList* _args=std::get_if<ValuesList>(&__args);
    if(not _args)
      errorEmitter("FuncCallNode args did not evaluate to a ValuesList but of type ",variantInnerTypeName(__args));
    
    const std::string fName=
      funcCallNode.fun;
    
    diagnostic("Going to call function ",fName,"\n");
    
    std::vector<std::shared_ptr<Value>> args;
    const Value& _fun=curEnv->at(fName);
    
    if(const HostFunction* fun=std::get_if<HostFunction>(&_fun))
      {
	args.reserve(_args->data.size());
	for(std::shared_ptr<Value> _arg : _args->data)
	  {
	    FuncArg& arg=unvariant<FuncArg>(*_arg);
	    
	    if(arg.name!="")
	      errorEmitter("trying to call a host function with named arg ",arg.name);
	    args.push_back(arg.value);
	  }
      }
    else
      if(const Function* fun=std::get_if<Function>(&_fun))
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
	     &fName](const bool insertNamed)
	    {
	      for(const std::shared_ptr<Value>& __arg : _args->data)
		{
		  FuncArg& _arg=unvariant<FuncArg>(*__arg);
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
			    std::get_if<ValuesList>(&*args.back())->data.push_back(arg);
			  else
			    errorEmitter("passing too many arguments to non-variadic function ",fName);
			}
		      else
			args[iPar]=arg;
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
	if(const StructDef* _sd=std::get_if<StructDef>(&_fun))
	  {
	    const StructDefNode& sd=*_sd->structDefNode;
	    const StructMemberListNode& members=sd.members;
	    
	    const size_t nMembers=members.list.size();
	    args.resize(nMembers);
	    
	    auto insertNamedOrUnnamed=
	      [&args,
	       &_args,
	       &members,
	       &fName](const bool insertNamed)
	      {
		for(const std::shared_ptr<Value>& __arg : _args->data)
		  {
		    FuncArg& _arg=unvariant<FuncArg>(*__arg);
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
			  args[iMember]=arg;
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
	  errorEmitter("trying to call a non-callable variable of kind: ",variantInnerTypeName(_fun));
    
    return callFunction(_fun,args);
  }
  
  EvalResult operator()(const FuncDefNode& funcDefNode)
  {
    const std::string& name=
      funcDefNode.name;
    
    if(curEnv->find(name))
      errorEmitter("Redefining a function which has a name ",name," already defined");
    
    (*curEnv)[name]=(*this)(*funcDefNode.fun).asRhs();
    
    return {};
  }
  
  EvalResult operator()(const FuncNode& funcNode)
  {
    return
      {.value=Function{.funcNode=&funcNode,
		       .env=curEnv}};
  }
  
  EvalResult operator()(const ASTNodesNode& astNodesNode)
  {
    EnvironmentWindUnwinder windUnwind=pushNewEnvironment();
    
    for(auto& n : astNodesNode.list)
      {
	if(std::get_if<ReturnNode>(&*n))
	  return std::visit(*this,*n);
	else
	  std::visit(*this,*n);
      }
    
    return {};
  }
  
  EvalResult operator()(const ReturnNode& returnNode)
  {
    return std::visit(*this,*returnNode.arg);
  }
  
  EvalResult operator()(const PrePostfixOpNode& prePostfixOpNode)
  {
    return{.value=prePostfixOpNode.op(std::visit(*this,*prePostfixOpNode.arg).asLhs())};
  }
  
  EvalResult operator()(const AssignNode& assignNode)
  {
    EvalResult _lhs=std::visit(*this,*assignNode.lhs);
    
    _lhs.set(assignNode.mode,std::visit(*this,*assignNode.rhs).asRhs());
    
    return _lhs;
  }
  
  EvalResult operator()(const UnOpNode& unOpNode)
  {
    return {.value=unOpNode.op(std::visit(*this,*unOpNode.arg).asRhs())};
  }
  
  EvalResult operator()(const BinOpNode& binOpNode)
  {
    return {.value=binOpNode.op(std::visit(*this,*binOpNode.arg1).asRhs(),std::visit(*this,*binOpNode.arg2).asRhs())};
  }
};

