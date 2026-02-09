// Nissa Interface for Correlator Evaluation

/**
  The logic

  Creation of the AST
  -------------------
  Within a script, one can define sources, operators, and functions,
  which are multiplied or summed. This are represented in terms of an
  AST, including function definitions and calls, assignement, etc.
  
  Execution of the AST
  -------------------
  When the AST is exectued, variables are
  
  Creation of the lines
  -------------------------------
   (contraction for now) triggers the
  creation of the lines
  
**/

#include <filesystem>
#include <memory>
#include <variant>

#include <nissa.hpp>

#include "nASTy.hpp"

using namespace nissa;

using CAction=
  std::tuple<std::string,std::function<void()>,std::vector<std::string>>;

/// Creates a source
struct SourcePars
{
  static inline int glbId;
  
  int id;
  
  rnd_t noiseType;
  
  bool spinDilute;
  
  bool colDilute;
  
  SourcePars(const rnd_t& noiseType=RND_Z4,
	     const bool& spinDilute=true,
	     const bool& colDilute=false) :
    id{glbId++},
    noiseType{noiseType},
    spinDilute{spinDilute},
    colDilute{colDilute}
  {
  }
};

/// Select a time
struct TimeSelectPars
{
  int t;
};

using ASTNodePars=
  std::variant<std::monostate,
		 SourcePars,
		 TimeSelectPars>;
/////////////////////////////////////////////////////////////////

enum OpType{SUM,PROD};

struct ASTNodeOp
{
  std::variant<ASTNodePars,std::pair<OpType,std::array<std::shared_ptr<ASTNodeOp>,2>>> subNodes;
  
  template <typename A>
  ASTNodeOp(const A& a) :
    subNodes(a)
  {
  }
};

std::shared_ptr<ASTNodeOp> source(const rnd_t& noiseType=RND_Z4,
				  const bool& spinDilute=true,
				  const bool& colDilute=false)
{
  return std::make_shared<ASTNodeOp>(SourcePars{noiseType,spinDilute,colDilute});
}

std::shared_ptr<ASTNodeOp> timeSelect(const int& t)
{
  return std::make_shared<ASTNodeOp>(TimeSelectPars{t});
}

// #define PROVIDE_AST_BINOP(A,B)						\
//   std::shared_ptr<ASTNodeOp> operator A(const std::shared_ptr<ASTNodeOp>& op1, \
// 					const std::shared_ptr<ASTNodeOp>& op2) \
//   {									\
//     return std::make_shared<ASTNodeOp>(std::make_pair(B,std::array<std::shared_ptr<ASTNodeOp>,2>{op1,op2})); \
//   }

// PROVIDE_AST_BINOP(+,SUM);
// PROVIDE_AST_BINOP(*,PROD);

// #undef PROVIDE_AST_BINOP

// enum LineT{QuarkLine,LepLine};

// struct Statement
// {
//   std::string out;
  
//   LineT lineT;
  
//   std::function<void(void)> eval;
  
//   std::vector<std::string> in;
// };

// Statement operator*(const TimeSelectPars& pars,
// 		    const Statement& oth)
// {
//   return {"DeltaT*"+oth.out,
//     oth.lineT,
//     [](){},
//     {oth.out}};
// }

int main()
{
  std::vector<char> ext;
  if(const char* path="/home/francesco/QCD/SORGENTI/nissa/test/projects/nice/input.cpp";
     std::filesystem::exists(path))
    {
      const size_t exs=
	std::filesystem::file_size(path);
      ext.resize(exs+1);
      
      if(FILE* file=fopen(path,"r"))
	{
	  if(const size_t n=fread(&ext[0],1,exs,file);n!=exs)
	    errorEmitter("expected ",exs," obtained ",n);
	  ext[exs]='\0';
	  fclose(file);
	}
      else
	errorEmitter("unable to read ",path);
    }
  else
    errorEmitter("file ",path," does not exists");
  
  pp::internal::verbose=true;
  
  const auto parser=
    getParser();
  
  const auto parseTree=
    createParseTree(parser,ext.data());
  
  const auto ptExecutor=
    getParseTreeExecuctor(parser.actionStrings);
  
  const auto ast=
    ptExecutor.execParseTree(parseTree);
  
  /////////////////////////////////////////////////////////////////
  
  Evaluator ev;
  
  pp::internal::diagnostic("Come on\n");
  
  std::visit(ev,*ast);
  
  // const auto eta=
  //   source();
  
  // MASTER_PRINTF("%d %d\n",SourcePars::glbId,std::get_if<SourcePars>(&eta->pars)->id);
  
  // auto p=
  //   timeSelect(0)*eta;
  
  // auto pr=
  //   [](const auto& pr,
  //      const auto& v)->void
  //   {
  //     auto npr=
  // 	Overload
  // 	{[](ASTNodePars pars)
  // 	{
  // 	  std::visit(Overload
  // 		     {[](const std::monostate& dum)
  // 		     {
  // 		       CRASH("monostate should never appear");
  // 		     },
  // 			 [](const TimeSelectPars& tPars)
  // 			 {
  // 			   MASTER_PRINTF("Selecting time: %d\n",tPars.t);
  // 			 },
  // 			 [](const SourcePars& sPars)
  // 			 {
  // 			   MASTER_PRINTF("Source with diluted spin: %d\n",sPars.spinDilute);
  // 			 }},pars);
  // 	},
  // 	 [&pr](const std::pair<OpType,std::array<std::shared_ptr<ASTNodeOp>,2>>& leafs)
  // 	 {
  // 	   MASTER_PRINTF("%s of: \n",(leafs.first==SUM)?"sum":"prod");
  // 	   pr(pr,leafs.second[0]->subNodes);
  // 	   pr(pr,leafs.second[1]->subNodes);
  // 	 }
  // 	};
      
  //     std::visit(npr,v);
  //   };
  
  // pr(pr,p->subNodes);
  
  return 0;
}
