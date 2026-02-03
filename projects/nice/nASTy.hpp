#include <cmath>

#include <parsePact.hpp>

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
    "             | \"return\" expression_statement [funcReturn($1)]"
    "             ;"
    "   functionDefinition : \"fun\" identifier \"\\(\" variadic_parameter_list \"\\)\" compound_statement [funcDef($1,$3,$5)]"
    "                      ;"
    "   variadic_parameter_list: [emptyParList]"
    "                          | parameterList [return($0)]"
    "                          | variadic_parameter [emptyVariadicParList($0)]"
    "                          | parameterList \",\" variadic_parameter [makeVariadicParList($0)]"
    "                          ;"
    "   variadic_parameter: \"\\.\\.\\.\" [return(1)]"
    "                     | \"&\" \"\\.\\.\\.\" [return(2)]"
    "                     ;"
    "   parameterList: parameter [firstParOfList]"
    "                 | parameterList \",\" parameter [appendParToList($0,$2)]"
    "                 ;"
    "   parameter: neededParameter [return($0)]"
    "            | neededParameter \"=\" expression [addParDefault($0,$2)]"
    "            ;"
    "   neededParameter: identifier [parCreate(0,$0)]"
    "                  | \"&\" identifier [parCreate(1,$1)]"
    "                  ;"
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
    "                       ;"
    "    funcArgsList : [emptyArgList]"
    "                 | nonemptyArgsList [return($0)]"
    "                 ;"
  "      nonemptyArgsList : arg [firstArgOfList]"
    "                     | funcArgsList \",\" arg [appendArgToList($0,$2)]"
    "                     ;"
    "    arg : \"\\.\" identifier \"=\" expression [createArg($1,$3)]"
    "        | expression [createArg(\"\",$0)]"
    "        ;"
    "    primary_expression : identifier [return($0)]"
    "                       | \"[0-9]+\" [convToInt($0)]"
    "                       | \"([0-9]+(\\.[0-9]*)?|(\\.[0-9]+))((e|E)(\\+|\\-)?[0-9]+)?\" [convToFloat($0)]"
    "                       | \"\\\"[^\\\"]*\\\"\" [convToStr($0)]"
    "                       | \"\\(\" expression \"\\)\" %precedence BRACKETS [return($1)]"
    "                       | \"lambda\" \"\\(\" variadic_parameter_list \"\\)\" compound_statement [lambdaFuncDef($2,$4)]"
    "                       ;"
    "    assign_expression : postfix_expression \"=\" expression %precedence \"=\" [unaryAssign($0,$2)]"
    "                      | postfix_expression \"\\*=\" expression [unaryProdAssign($0,$2)]"
    "                      | postfix_expression \"/=\" expression [unaryDivAssign($0,$2)]"
    "                      | postfix_expression \"\\+=\" expression [unarySumAssign($0,$2)]"
    "                      | postfix_expression \"\\-=\" expression [unaryDiffAssign($0,$2)]"
    "                      ;"
    "    identifier : \"[a-zA-Z_][a-zA-Z0-9_]*\" [convToId($0)]"
    "               ;"
    "}";
  
  return
    pp::createGrammar(cGrammar);
}

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

struct BinOpNode;

struct IdNode;

struct AssignNode;

struct ReturnNode;

struct FuncParNode;

struct FuncArgNode;

struct FuncArgListNode;

struct FuncParListNode;

struct SubscribeNode;

using ASTNode=
  std::variant<ASTNodesNode,
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
	       ReturnNode,
	       ValueNode,
	       AssignNode>;

struct ValueRef;

struct Function;

struct HostFunction;

struct ValuesList;

using Value=
  std::variant<std::monostate,std::string,int,double,Function,HostFunction,ValueRef,ValuesList>;

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
  const FuncNode* fun;
  
  Environment* env;
};

struct HostFunction :
  std::function<Value(std::vector<std::shared_ptr<Value>>&&)>
{
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

struct FuncArgListNode
{
  std::vector<FuncArgNode> list;
};

struct SubscribeNode
{
  std::shared_ptr<ASTNode> base;
  
  std::shared_ptr<ASTNode> subscr;
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
};

struct FuncParListNode
{
  std::vector<FuncParNode> list;
  
  auto iParOfName(const std::string& name) const
  {
    auto it=
      std::find_if(list.begin(),
		   list.end(),
		   [&name](const FuncParNode& par)
		   {
		     return par.name==name;
		   });
    
    return std::distance(list.begin(),it);
  }
  
  enum VariadicMode{NONE,BY_VALUE,BY_REF};
  
  VariadicMode variadicMode;
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

struct AssignNode
{
  std::shared_ptr<ASTNode> lhs;
  
  std::shared_ptr<ASTNode> rhs;
  
  std::function<Value(Value&,const Value&)> op;
};

struct PrePostfixOpNode
{
  std::shared_ptr<ASTNode> arg;
  
  std::function<Value(Value&)> op;
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
  Environment* parent;
  
  std::unordered_map<std::string,std::shared_ptr<Value>> varTable;
  
  Environment* uppermost()
  {
    Environment* res=this;
    
    while(Environment* par=res->parent)
      res=par;
    
    return res;
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
      pp::internal::errorEmitter("Trying to use uninitialized variable ",name);
    
    return *f;
  }
  
  std::shared_ptr<Value> getRefOrInsert(const std::string& name)
  {
    if(auto f=find(name))
      return f;
    else
      return varTable[name]=std::make_shared<Value>(std::monostate{});
  }
  
  Value& operator[](const std::string& name)
  {
    return *getRefOrInsert(name);
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
  
  Environment(Environment* parent=nullptr) :
    parent(parent)
  {
  }
};

template <typename T>
T& fetch(std::shared_ptr<ASTNode>& subNode,
	 const char* comm=nullptr)
{
  using namespace pp::internal;
  
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
  //  PROVIDE_ACTION_WITH_N_SYMBOLS("firstStatement",1,return std::make_shared<ASTNode>(ASTNodesNode{.list{subNodes[0]}}));
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
  PROVIDE_ACTION_WITH_N_SYMBOLS("createArg",2,
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
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("parCreate",2,return std::make_shared<ASTNode>(FuncParNode{.name=fetch<IdNode>(subNodes,1).name,
											   .isRef=(bool)unvariant<int>(fetch<ValueNode>(subNodes,0).value)}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("addParDefault",2,const FuncParNode& in=fetch<FuncParNode>(subNodes,0);
				return std::make_shared<ASTNode>(FuncParNode{.name=in.name,.isRef=in.isRef,.def=subNodes[1]}));
#define PROVIDE_FUNC_LIST_ACTIONS(NAME)					\
  PROVIDE_ACTION_WITH_N_SYMBOLS("empty" #NAME "List",0,return std::make_shared<ASTNode>(Func ## NAME ## ListNode{})); \
  PROVIDE_ACTION_WITH_N_SYMBOLS("first" #NAME "OfList",1,return makeFirstOfList<Func ## NAME ## Node,Func ## NAME ## ListNode>(subNodes[0])); \
  PROVIDE_ACTION_WITH_N_SYMBOLS("append" #NAME "ToList",2,return appendToList<Func ##NAME ## Node,Func ## NAME ## ListNode>(subNodes[0],subNodes[1])); \
  //  PROVIDE_ACTION_WITH_N_SYMBOLS("merge" #NAME "Lists",2,return mergeLists<Func ## NAME ## ListNode>(subNodes[0],subNodes[1]))
  
  PROVIDE_FUNC_LIST_ACTIONS(Par);
  PROVIDE_FUNC_LIST_ACTIONS(Arg);
  
#undef PROVIDE_FUNC_LIST_ACTIONS
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("emptyVariadicParList",1,return std::make_shared<ASTNode>(FuncParListNode{.variadicMode=FuncParListNode::VariadicMode(unvariant<int>(fetch<ValueNode>(subNodes,0).value))}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("makeVariadicParList",2,
				fetch<FuncParListNode>(subNodes,0).variadicMode=FuncParListNode::VariadicMode(unvariant<int>(fetch<ValueNode>(subNodes,1).value));
				return subNodes[0]);
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("funcDef",3,
				return std::make_shared<ASTNode>(FuncDefNode{.name=fetch<IdNode>(subNodes,0).name,
									     .fun=std::make_shared<FuncNode>(fetch<FuncParListNode>(subNodes,1),
													     subNodes[2])}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("lambdaFuncDef",2,
				return std::make_shared<ASTNode>(FuncNode{fetch<FuncParListNode>(subNodes,0),
									  subNodes[1]}));
  
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("funcCall",2,return std::make_shared<ASTNode>(FuncCallNode{.fun=fetch<IdNode>(subNodes,0).name,
											   .args=fetch<FuncArgListNode>(subNodes,1)}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("funcReturn",1,return std::make_shared<ASTNode>(ReturnNode{.arg=subNodes[0]}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("subscribe",2,return std::make_shared<ASTNode>(SubscribeNode{.base=subNodes[0],.subscr=subNodes[1]}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("forStatement",4,return std::make_shared<ASTNode>(ForNode{.subNodes{subNodes}}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("ifStatement",2,return std::make_shared<ASTNode>(IfNode{.subNodes{subNodes}}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("ifElseStatement",3,return std::make_shared<ASTNode>(IfNode{.subNodes{subNodes}}));
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("unaryAssign",2,return std::make_shared<ASTNode>(AssignNode{.lhs=subNodes[0],
	  .rhs=subNodes[1],
	  .op=[](Value& lhs,
		 const Value& rhs)
	  {
	    return lhs=rhs;
	  }}));

#define PROVIDE_ASSIGN(NAME,SYMBOL)					\
  PROVIDE_ACTION_WITH_N_SYMBOLS("unary" #NAME "Assign",2,return std::make_shared<ASTNode>(AssignNode{.lhs=subNodes[0], \
	  .rhs=subNodes[1],						\
	  .op=[](Value& lhs,						\
		 const Value& rhs)					\
	  {								\
	    std::visit([](auto& lhs,					\
			  const auto& rhs)				\
	    {								\
	      if constexpr(requires {lhs SYMBOL ## =rhs;})		\
		lhs SYMBOL ## =rhs;					\
	      else							\
		pp::internal::errorEmitter("Cannot " #SYMBOL" the types:",typeid(lhs).name()," and ",typeid(rhs).name()); \
	    },lhs,rhs);							\
	    								\
	    return lhs;							\
	  }}))
  
  PROVIDE_ASSIGN(Prod,*);
  PROVIDE_ASSIGN(Div,/);
  PROVIDE_ASSIGN(Sum,+);
  PROVIDE_ASSIGN(Diff,-);
  
#undef PROVIDE_ASSIGN
  
#undef PROVIDE_ACTION_WITH_N_SYMBOLS
  
  return ParseTreeExecutor(providedActions,requiredActions);
}

struct Evaluator
{
  Environment env;
  
  int iLhs{};
  
  void prepareTopEnv()
  {
    env["M_PI"]=M_PI;
    
    env["print"]=
      HostFunction{
      [](std::vector<std::shared_ptr<Value>>&& args)->Value
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
      [](std::vector<std::shared_ptr<Value>>&& args)->Value
      {
	return ValuesList{args};
      }};
    
#define REGISTER_ARGLESS_HOST_FUNCTION(NAME)	\
    env[#NAME]=					\
      HostFunction{				\
      [](std::vector<std::shared_ptr<Value>>&& args)->Value	\
      {						\
	return NAME();				\
      }}
    
#define REGISTER_HOST_FUNCTION(NAME,ARGS...)				\
    env[#NAME]=								\
      HostFunction{							\
      [](std::vector<std::shared_ptr<Value>>&& args)->Value		\
      {									\
	using pp::internal::errorEmitter;				\
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
	return std::visit([](auto&&...args) ->Value			\
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
  
  Value maybeEvalAsLhs(std::shared_ptr<ASTNode> arg)
  {
    iLhs++;
    Value _lhs=std::visit(*this,*arg);
    iLhs--;
    
    return _lhs;
  }
  
  Value& evalAsLhs(std::shared_ptr<ASTNode> arg)
  {
    Value _lhs=maybeEvalAsLhs(arg);
    
    ValueRef* lhs=std::get_if<ValueRef>(&_lhs);
    
    if(not lhs)
      pp::internal::errorEmitter("arg does not eval to a l-value ref");
    
    return *lhs->ref;
  }
  
  Value operator()(const ValueNode& valueNode)
  {
    return valueNode.value;
  }
  
  Value operator()(const IdNode& idNode)
  {
    if(iLhs)
      return ValueRef{env.getRefOrInsert(idNode.name)};
    else
      return env.at(idNode.name);
  }
  
  Value operator()(const SubscribeNode& subscribeNode)
  {
    Value b=maybeEvalAsLhs(subscribeNode.base);
    
    ValueRef* vr=std::get_if<ValueRef>(&b);
    Value& v=vr?(*vr->ref):b;
    
    Value s=std::visit(*this,*subscribeNode.subscr);
    
    using pp::internal::errorEmitter;
    using pp::internal::variantInnerTypeName;
    
    if(ValuesList* vl=std::get_if<ValuesList>(&v))
      if(int* mi=std::get_if<int>(&s))
	{
	  int i=*mi;
	  
	  if(i<0)
	    errorEmitter("Trying to subscribe with negative index ",i);
	  
	  if(i>=vl->data.size())
	    errorEmitter("Trying to subscribe withindex ",i," greater than list size");
	  
	  std::shared_ptr<Value> r=vl->data[i];
	  
	  if(iLhs)
	    return (ValueRef){r};
	  else
	    return *r;
	}
      else
	errorEmitter("index of subscrition is not an integer but is of type: ",variantInnerTypeName(s));
    else
      errorEmitter("trying to subscribe non-list of type: ",variantInnerTypeName(b));
    
    return {};
  }
  
  Value operator()(const ForNode& forNode)
  {
    Evaluator subev{&env};
    
    for(std::visit(*this,*forNode.subNodes[0]);
    	std::visit([](const auto& v)
	{
	  if constexpr(std::is_convertible_v<decltype(v),bool>)
	    return (bool)v;
	  else
	    pp::internal::errorEmitter("Cannot convert the type to bool");
	  
	  return false;
	},std::visit(subev,*forNode.subNodes[1]));
	std::visit(subev,*forNode.subNodes[2]))
      std::visit(subev,*forNode.subNodes[3]);
    
    return std::monostate{};
  }
  
  Value operator()(const IfNode& ifNode)
  {
    Evaluator subev{&env};
    
    if(std::visit([](const auto& v)
	{
	  if constexpr(std::is_convertible_v<decltype(v),bool>)
	    return (bool)v;
	  else
	    pp::internal::errorEmitter("Cannot convert the type to bool");
	  
	  return false;
	},std::visit(subev,*ifNode.subNodes[0])))
      std::visit(subev,*ifNode.subNodes[1]);
    else
      if(ifNode.subNodes.size()>=2)
	std::visit(subev,*ifNode.subNodes[2]);
    
    return std::monostate{};
  }
  
  Value operator()(const FuncParNode&)
  {
    return std::monostate{};
  }
  
  Value operator()(const FuncParListNode&)
  {
    return std::monostate{};
  }
  
  Value operator()(const FuncArgNode&)
  {
    return std::monostate{};
  }
  
  Value operator()(const FuncArgListNode&)
  {
    return std::monostate{};
  }
  
  Value operator()(const FuncCallNode& funcCallNode)
  {
    using namespace pp::internal;
    
    const std::string name=
      funcCallNode.fun;
    
    diagnostic("Going to call function ",name.c_str(),"\n");
    
    const std::shared_ptr<Value> fv=env.find(name);
    if(not fv)
      errorEmitter("unable to find function ",name);
    
    if(const Function* _f=std::get_if<Function>(&*fv))
      {
	const Function& f=*_f;
	
	const FuncNode& fn=*f.fun;
	const FuncParListNode& pars=fn.pars;
	
	Environment* fe=f.env;
	if(f.env==nullptr)
	  fe=this->env.uppermost();
	
	Evaluator subev{fe};
	
	auto getValOrRef=
	  [this](const std::string& name,
		 const bool& isRef,
		 const std::shared_ptr<ASTNode>& expr)
	  {
	    std::shared_ptr<Value> v=
	      std::make_shared<Value>(std::visit(*this,*expr));
	    
	    if(isRef)
	      {
		if(ValueRef* vr=std::get_if<ValueRef>(&*v))
		  v=vr->ref;
		else
		  errorEmitter("argument ",name," has to be taken by ref but does not evalue to a ref");
	      }
	    
	    return v;
	  };
	
	auto insertInSubev=
	  [&subev,
	   &getValOrRef](const std::string& name,
			 const bool& isRef,
			 const std::shared_ptr<ASTNode>& expr)
	  {
	    auto r=
	      subev.env.varTable.try_emplace(name,getValOrRef(name,isRef,expr));
	    
	    if(not r.second)
	      errorEmitter("argument ",name," already specified");
	  };
	
	size_t iNextPositionalArg{};
	bool acceptingOnlyNamedArgs{};
	std::shared_ptr<Value> vaArgsHandle;
	const bool nFixedArgs=fn.pars.list.size();
	for(const FuncArgNode& ap : funcCallNode.args.list)
	  {
	    const bool isNamedArg=ap.name!="";
	    acceptingOnlyNamedArgs|=isNamedArg;
	    
	    if(isNamedArg)
	      {
		const std::string& name=
		  ap.name;
		
		const size_t iPar=
		  pars.iParOfName(name);
		if(iPar==pars.list.size())
		  errorEmitter("trying to pass argument ",name," not expected by the function");
		
		const FuncParNode& parDef=
		  pars.list[iPar];
		
		printf("Plugging named par: %s\n",name.c_str());
		insertInSubev(name,parDef.isRef,ap.expr);
	      }
	    else
	      {
		if(acceptingOnlyNamedArgs)
		  errorEmitter("Positional arg passed after named ones");
		
		if(iNextPositionalArg>=nFixedArgs)
		  {
		    if(not vaArgsHandle)
		      {
			vaArgsHandle=std::make_shared<Value>(std::in_place_type<ValuesList>);
			
			if(pars.variadicMode!=FuncParListNode::VariadicMode::NONE)
			  errorEmitter("Trying to pass ",iNextPositionalArg," arguments, more than those expected by a non-variadic function");
		      }
		    
		    ValuesList& vaArgs=
		      std::get<ValuesList>(*vaArgsHandle);
		    
		    const bool isRef=
		      pars.variadicMode==FuncParListNode::VariadicMode::BY_REF;
		    
		    vaArgs.data.emplace_back(getValOrRef("Variadic argument number "+std::to_string(iNextPositionalArg-nFixedArgs),isRef,ap.expr));
		  }
		else
		  {
		    const FuncParNode& par=
		      pars.list[iNextPositionalArg];
		    printf("Plugging pos par: %s\n",par.name.c_str());
		    insertInSubev(par.name,par.isRef,ap.expr);
		  }
	      }
	    
	    iNextPositionalArg++;
	  }
	
	if(vaArgsHandle)
	  if(not subev.env.varTable.try_emplace("__VA_ARGS__",vaArgsHandle).second)
	    errorEmitter("__VA_ARGS__ already defined as a parameter");
	
	// diagnostic("Calling function, specified arguments:\n");
	// subev.env.print();
	
	// Put possible default pars
	for(auto& [aName,isRef,optDef] : fn.pars.list)
	  if(not subev.env.varTable.contains(aName))
	    if(optDef)
	      subev.env[aName]=std::visit(*this,*optDef);
	    else
	      errorEmitter("parameter \"",aName,"\" with no default value unspecified when calling the function");
	  else
	    {}
	
	return std::visit(subev,*fn.body);
      }
    else
      if(const HostFunction* hf=std::get_if<HostFunction>(&*fv))
	{
	  std::vector<std::shared_ptr<Value>> evArgs;
	  evArgs.reserve(funcCallNode.args.list.size());
	  for(const FuncArgNode& fan : funcCallNode.args.list)
	    evArgs.push_back(std::make_shared<Value>(std::visit(*this,*fan.expr)));
	  
	  return (*hf)(std::move(evArgs));
	}
      else
	errorEmitter("Variable is of type ",variantInnerTypeName(*fv)," not a ",typeid(Function).name());
    
    return std::monostate{};
  }
  
  Value operator()(const FuncDefNode& funcDefNode)
  {
    const std::string& name=
      funcDefNode.name;
    
    if(env.find(name))
      pp::internal::errorEmitter("Redefining a function which has a name ",name," already defined");
    
    env[name]=(*this)(*funcDefNode.fun);
    
    return std::monostate{};
  }
  
  Value operator()(const FuncNode& funcNode)
  {
    return
      Function{.fun=&funcNode,
	       .env=&env};
  }
  
  Value operator()(const ASTNodesNode& astNodesNode)
  {
    /// Creates a subev only if there is already one above
    std::shared_ptr<Evaluator> _subev;
    Evaluator *subev=this;
    if(env.parent!=nullptr)
      {
	_subev=std::make_shared<Evaluator>(&env);
	subev=_subev.get();
      }
    
    for(auto& n : astNodesNode.list)
      {
	if(std::get_if<ReturnNode>(&*n))
	  return std::visit(*subev,*n);
	else
	  std::visit(*subev,*n);
      }
    
    return {};
  }
  
  Value operator()(const ReturnNode& returnNode)
  {
    return std::visit(*this,*returnNode.arg);
  }
  
  Value operator()(const PrePostfixOpNode& prePostfixOpNode)
  {
    return
      prePostfixOpNode.op(evalAsLhs(prePostfixOpNode.arg));
  }
  
  Value operator()(const AssignNode& assignNode)
  {
    return
      assignNode.op(evalAsLhs(assignNode.lhs),std::visit(*this,*assignNode.rhs));
  }
  
  Value operator()(const UnOpNode& unOpNode)
  {
    return unOpNode.op(std::visit(*this,*unOpNode.arg));
  }
  
  Value operator()(const BinOpNode& binOpNode)
  {
    return binOpNode.op(std::visit(*this,*binOpNode.arg1),std::visit(*this,*binOpNode.arg2));
  }
};

