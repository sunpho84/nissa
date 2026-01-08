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
    "   statement : expression_statement [return]"
    "             | compound_statement [return]"
    "             | forStatement [return]"
    "             | ifStatement [return]"
    "             | functionDefinition [return]"
    "             | \"return\" expression_statement [funcReturn(1)]"
    "             ;"
    "   functionDefinition : \"fun\" identifier \"\\(\" expressionList \"\\)\" compound_statement [funcDef(1,3,5)]"
    "                      ;"
    "   forStatement: \"for\" \"\\(\" forInit \";\" forCheck \";\" forIncr \"\\)\" statement [forStatement(2,4,6,8)]"
    "               ;"
    "   forInit: expression [return]"
    "          |"
    "          ;"
    "   forCheck: expression [return]"
    "           |"
    "           ;"
    "   forIncr: expression [return]"
    "          |"
    "          ;"
    "   ifStatement: \"if\" \"\\(\" expression \"\\)\" statement %precedence lowerThanElse [ifStatement(2,4)]"
    "              | \"if\" \"\\(\" expression \"\\)\" statement \"else\" statement [ifElseStatement(2,4,6)]"
    "              ;"
    "    compound_statement : \"{\" statements \"}\" [return(1)]"
    "                       ;"
    "    statements : [createStatements]"
    "               | statements statement [appendStatement]"
    "               ;"
    "    expression_statement : expression \";\" [return(0)]"
    "                         ;"
    "    expression : postfix_expression %precedence \"\\+\\+\" [return(0)]"
    "               | expression \"\\+\" expression [binarySum(0,2)]"
    "               | expression \"\\-\" expression [binaryDiff(0,2)]"
    "               | expression \"\\*\" expression [binaryProd(0,2)]"
    "               | expression \"/\" expression [binaryDiv(0,2)]"
    "               | expression \"%\" expression [binaryMod(0,2)]"
    "               | \"\\+\" expression %precedence UNARY_ARITHMETIC [unaryPlus(1)]"
    "               | \"\\-\" expression %precedence UNARY_ARITHMETIC [unaryMinus(1)]"
    "               | \"!\" expression [unaryNot(1)]"
    "               | expression \"<\" expression [binarySmaller(0,2)]"
    "               | expression \">\" expression [binaryGreater(0,2)]"
    "               | expression \"<=\" expression [binarySmallerEqual(0,2)]"
    "               | expression \">=\" expression [binaryGreaterEqual(0,2)]"
    "               | expression \"==\" expression [binaryCompare(0,2)]"
    "               | expression \"!=\" expression [binaryInequal(0,2)]"
    "               | expression \"&&|and\" expression [binaryAnd(0,2)]"
    "               | expression \"\\|\\||or\" expression [binaryOr(0,2)]"
    "               | assign_expression [return]"
    "               | identifier \"\\*=\" expression [unaryProdAssign]"
    "               | identifier \"/=\" expression [unaryDivAssign]"
    "               | identifier \"\\+=\" expression [unarySumAssign]"
    "               | identifier \"\\-=\" expression [unaryDiffAssign]"
    "               ;"
    "    postfix_expression : primary_expression [return(0)]"
    "                       | postfix_expression \"\\+\\+\" [unaryPostfixIncrement(0)]"
    "                       | postfix_expression \"\\-\\-\" [unaryPostfixDecrement(0)]"
    // "                       | postfix_expression \"\\(\" \"\\)\" %precedence FUNCTION_CALL [emptyFuncCall(0)] "
    "                       | postfix_expression \"\\(\" funcArgsList \"\\)\" %precedence FUNCTION_CALL [funcCall(0,2)] "
    "                       | postfix_expression \"\\[\" expression \"\\]\" [subscribe(0,2)] "
    "                       ;"
    "    funcArgsList : expressionList [return]"
    "                 | namedArgList [return]"
    "                 | expressionList \",\" namedArgList [immediateSum(0,2)]"
    "                 ;"
    "    namedArgList : namedArg [firstExprOfList]"
    "                 | namedArgList \",\" namedArg [appendStatement(0,2)]"
    "                 ;"
    "    namedArg : \"\\.\" identifier \"=\" expression [namedArg(1,3)]"
    "             ;"
    "    primary_expression : identifier [return]"
    "                       | \"[0-9]+\" [convToInt]"
    "                       | \"([0-9]+(\\.[0-9]*)?|(\\.[0-9]+))((e|E)(\\+|\\-)?[0-9]+)?\" [convToFloat]"
    "                       | \"\\\"[^\\\"]*\\\"\" [convToStr]"
    "                       | \"\\(\" expression \"\\)\" %precedence BRACKETS [return(1)]"
    "                       ;"
    "    assign_expression : postfix_expression \"=\" expression %precedence \"=\" [unaryAssign(0,2)]"
    "                      ;"
    // "    lhs : identifier [getPtrOfId]"
    // "        | lhs \"[\" expression \"]\" [subscribe(0,2)]"
    // "        ;"
    "    expressionList : [createStatements]"
    "                   | nonEmptyExpressionList %precedence \",\" [return]"
    "                   ;"
    "    nonEmptyExpressionList : expression [firstExprOfList]"
    "                           | nonEmptyExpressionList \",\" expression [appendStatement(0,2)]"
    "                           ;"
    "    identifier : \"[a-zA-Z_][a-zA-Z0-9_]*\" [convToId]"
    "               ;"
    "}";
  
  return
    pp::createGrammar(cGrammar);
}

struct Environment;

struct ASTNodesNode;

struct ValueNode;

struct UnOpNode;

struct BinOpNode;

struct IdNode;

struct AssignNode;

struct ReturnNode;

using ASTNode=
  std::variant<ASTNodesNode,
	       // UnOpNode,
	       // BinOpNode,
	       // ForNode,
	       // IfNode,
	       IdNode,
	       // UPlusNode,
	       // UMinusNode,
	       // UNotNode,
	       // UPostfixIncrementNode,
	       // UPostfixDecrementNode,
	       // SumNode,
	       // DiffNode,
	       // ProdNode,
	       // DivNode,
	       // ModNode,
	       // SmallerNode,
	       // GreaterNode,
	       // SmallerEqualNode,
	       // GreaterEqualNode,
	       // CompareNode,
	       // OrNode,
	       // AndNode,
	       // InequalNode,
	       // FuncDefNode,
	       // FuncCallNode,
	       // SubscribeNode,
	       // NamedArgNode,
	       ReturnNode,
	       ValueNode,
	       AssignNode
  >;

struct ValueRef;

using Value=
  std::variant<std::monostate,std::string,int,ValueRef>;

struct ValueRef
{
  std::shared_ptr<Value> ref;
};

struct IdNode
{
  std::string name;
};

struct ValueNode
{
  Value value;
};

struct ASTNodesNode
{
  std::vector<std::shared_ptr<ASTNode>> subNodes;
};

struct AssignNode
{
  std::shared_ptr<ASTNode> lhs;
  
  std::shared_ptr<ASTNode> rhs;
  
  std::function<Value(Value&,const Value&)> op;
};

struct UnOpNode
{
  std::shared_ptr<ASTNode> arg;
};

struct BinOpNode
{
  std::shared_ptr<ASTNode> arg1;
  
  std::shared_ptr<ASTNode> arg2;
};

struct ReturnNode
{
  std::shared_ptr<ASTNode> arg;
};

struct Environment
{
  Environment* parent;
  
  std::unordered_map<std::string,std::shared_ptr<Value>> varTable;
  
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
      return varTable[name]=std::make_shared<Value>();
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
T& fetch(std::vector<std::shared_ptr<ASTNode>>& subNodes,
	 const size_t& i)
{
  using namespace pp::internal;
  
  if(const size_t n=subNodes.size();n<i)
    errorEmitter(n," nodes received, aksed for node #",i);
  
  T* s=
    std::get_if<T>(&*subNodes[i]);
  
  if(not s)
    errorEmitter("subNode ",i," is not of the required type ",typeid(T).name()," but is of type ",variantInnerTypeName(*subNodes[i]));
  
  return *s;
}

inline auto getParseTreeExecuctor()
{
  using namespace pp::internal;
  
  pp::ParseTreeExecutor<ASTNode,ValueNode> ptExecutor;
  
#define ENSURE_N_SYMBOLS(N)			\
  if(subNodes.size()!=N)			\
    errorEmitter("expecting " #N " symbols, obtained ",subNodes.size())
  
#define PROVIDE_ACTION_WITH_N_SYMBOLS(NAME,				\
				      N,				\
				      BODY...)				\
  ptExecutor.actions[NAME]=						\
    [](std::vector<std::shared_ptr<ASTNode>>& subNodes) -> std::shared_ptr<ASTNode> \
    {									\
      ENSURE_N_SYMBOLS(N);						\
      									\
      BODY;								\
    }
  
  PROVIDE_ACTION_WITH_N_SYMBOLS("createStatements",0,return std::make_shared<ASTNode>(ASTNodesNode{}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("appendStatement",2,
				if(ASTNodesNode* l=std::get_if<ASTNodesNode>(&*(subNodes[0]));l==nullptr)
				  errorEmitter("first argument is not a list of statement: ",std::visit([](const auto& x){return typeid(decltype(x)).name();},*(subNodes[0])));
				else
				  l->subNodes.emplace_back(subNodes[1]);
				return subNodes[0]);
  PROVIDE_ACTION_WITH_N_SYMBOLS("convToInt",1,return std::make_shared<ASTNode>(ValueNode{atoi(unvariant<std::string>(fetch<ValueNode>(subNodes,0).value).c_str())}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("return",1,return subNodes[0]);
  PROVIDE_ACTION_WITH_N_SYMBOLS("convToId",1,
				return std::make_shared<ASTNode>(IdNode{.name=unvariant<std::string>(fetch<ValueNode>(subNodes,0).value)}));
  PROVIDE_ACTION_WITH_N_SYMBOLS("unaryAssign",2,return std::make_shared<ASTNode>(AssignNode{.lhs=subNodes[0],
											    .rhs=subNodes[1],
											    .op=[](Value& lhs,
												   const Value& rhs)
											    {
											      return lhs=rhs;
											    }}));
  
  return ptExecutor;
}

struct Evaluator
{
  Environment env;
  
  bool isLhs=false;
  
  Value operator()(const ValueNode& valueNode)
  {
    return valueNode.value;
  }
  
  Value operator()(const IdNode& idNode)
  {
    printf("Going to evaluate idNode %s as lhs: %d\n",idNode.name.c_str(),isLhs);
    
    if(isLhs)
      return ValueRef{env.getRefOrInsert(idNode.name)};
    else
      return env.at(idNode.name);
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
    
    for(auto& n : astNodesNode.subNodes)
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
  
  Value operator()(const AssignNode& assignNode)
  {
    Value _lhs=std::visit(Evaluator(&env,true),*assignNode.lhs);
    
    ValueRef* lhs=std::get_if<ValueRef>(&_lhs);
    
    if(not lhs)
      pp::internal::errorEmitter("Lhs does not eval to a l-value ref");
    
    return
      assignNode.op(*lhs->ref,std::visit(*this,*assignNode.rhs));
  }
};

