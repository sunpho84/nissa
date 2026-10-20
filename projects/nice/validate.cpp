#include "corrParts.hpp"
#include "correlations.hpp"
#include "correlationTaskGraph.hpp"
#include "graphCompiler.hpp"
#include "schedulingStrategies.hpp"
#include "legacyInputWriter.hpp"
#include "niceSetup.hpp"
#include "nice.hpp"
#include "nASTy.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace nissa::nice;
bool containsSelection(const Line& line,
                       const int time)
{
  const auto* extender=
    std::get_if<Extender>(&line.data);
  
  if(not extender)
    return false;
  
  if(extender->op.tSelect==time and std::get_if<Gamma>(&extender->op.pars))
    return true;
  
  for(const auto& [weight,input] : extender->rhs)
    {
      (void)weight;
      if(containsSelection(input,time))
        return true;
    }
  
  return false;
}

void validateCorrParts()
{
  const Source eta(0);
  const Source copy=
    eta;
  const Source independent(0);
  assert(copy==eta);
  assert(independent!=eta);
  
  const Prop prop{.kappa=0.1394267,.mass=0.00066690,.r=-1,
                  .residue=1e-20};
  const Prop reversed=
    prop.dag();
  assert(reversed.r==1);
  assert(reversed.dag()==prop);
  
  const Momentum momentum{1,2,3};
  assert(momentum[0]==1 and momentum[2]==3);
  assert(((2.0*momentum).asArray()==std::array<double,3>{2,4,6}));
  assert(Phase{.mom=momentum}.dag().mom==Momentum({-1,-2,-3}));
  
  const Oper op=
    prop*Phase{.mom=momentum}*Gamma{.iGamma=5};
  const Oper reversedOp=
    op.dag();
  assert(reversedOp.dag()==op);
  assert((2.0*op).weight==2.0);
  
  const Line original=
    op*DeltaT{.t=5}*eta;
  const Line named(original,"first");
  const Line renamed(original,"second");
  assert(named==renamed);
  assert(named.describe()==renamed.describe());
  assert(containsSelection(original,5));
  
  const Line sum=
    original+2.0*original;
  const auto* extender=
    std::get_if<Extender>(&sum.data);
  assert(extender and extender->rhs.size()==2);
  assert(extender->rhs[0].first==1.0);
  assert(extender->rhs[1].first==2.0);
  assert(addLines(original,original,false)!=original);
}

void validateTaskGraph(int argc,char** argv)
{
  const Source source(0);
  CorrelationTaskGraph graph({Task{.id=0,.operation=source,.inputs{},.labels{"source"}},
                   Task{.id=1,
                        .operation=ExtendSelectingTPars{Pars(Gamma{.iGamma=1})},
                        .inputs{{0,1.0}},.labels{"line"}}},
                  {OutputBlock{.name="empty",.gammaList{},.traces{}}});
  graph.validate();
  assert(graph.tasks.size()==2);
  assert(graph.outputs.front().name=="empty");
  
  if(argc>1 and std::string(argv[1])=="invalid-cycle")
    {
      graph.tasks[1].inputs[0].task=1;
      graph.validate();
    }
  
  if(argc>1 and std::string(argv[1])=="invalid-graph-output")
    {
      graph.outputs.front().traces.push_back(5);
      graph.validate();
    }
}

void validateGraphCompiler(int argc,char** argv)
{
  const Source source(0);
  const Source copied=source;
  const Source independent(0);
  const Line backward=Prop{.r=1}*source;
  const Line forward=Gamma{.iGamma=1}*source;
  const Line forwardAlias(forward,"sameForward");
  const Line other=Gamma{.iGamma=1}*independent;
  const Line weighted(Pars(Gamma{.iGamma=0}),
                      {{2.0,forward},{3.0,Gamma{.iGamma=1}*copied}});
  const Line selected(Pars(DeltaT{.t=4}),{{1.0,Line(source)}});
  
  Correlations requests;
  Correlation& first=requests.addBlock("first");
  first.addGammas(5,5);
  first.addTrace(backward,forward);
  first.addTrace(weighted,selected);
  Correlation& second=requests.addBlock("second");
  second.addGammas(5,5);
  second.addTrace(backward,forwardAlias);
  second.addTrace(forward,backward);
  second.addTrace(backward,other);
  Correlation& empty=requests.addBlock("empty");
  empty.addGammas(0,0);
  
  const CorrelationTaskGraph graph=buildTaskGraph(requests);
  const CorrelationTaskGraph again=buildTaskGraph(requests);
  graph.validate();
  assert(graph.tasks==again.tasks);
  assert(graph.outputs==again.outputs);
  assert(graph.outputs.size()==3);
  assert(graph.outputs[2].name=="empty" and graph.outputs[2].traces.empty());
  
  const TaskId sharedTrace=graph.outputs[0].traces[0];
  assert(sharedTrace==graph.outputs[1].traces[0]);
  assert(sharedTrace!=graph.outputs[1].traces[1]);
  assert(graph.outputs[1].traces[1]!=graph.outputs[1].traces[2]);
  const Task& trace=graph.tasks[sharedTrace];
  assert(std::holds_alternative<TraceContraction>(trace.operation));
  assert(trace.inputs.size()==2);
  assert(trace.inputs[0].task!=trace.inputs[1].task);
  assert(trace.inputs[0].weight==1.0 and trace.inputs[1].weight==1.0);
  assert(std::find(graph.tasks[trace.inputs[1].task].labels.begin(),
                   graph.tasks[trace.inputs[1].task].labels.end(),
                   "sameForward")!=graph.tasks[trace.inputs[1].task].labels.end());
  
  const Task& weightedTrace=graph.tasks[graph.outputs[0].traces[1]];
  const Task& weightedTask=graph.tasks[weightedTrace.inputs[0].task];
  assert(weightedTask.inputs.size()==1);
  assert(weightedTask.inputs[0].task==trace.inputs[1].task);
  assert(weightedTask.inputs[0].weight==5.0);
  const Task& selectedTask=graph.tasks[weightedTrace.inputs[1].task];
  assert(std::get<ExtendSelectingTPars>(selectedTask.operation).tSelect==4);
  assert(selectedTask.inputs[0].task==graph.tasks[trace.inputs[0].task].inputs[0].task);
  
  const Task& independentTrace=graph.tasks[graph.outputs[1].traces[2]];
  assert(independentTrace.inputs[1].task!=trace.inputs[1].task);
  
  if(argc>1 and std::string(argv[1])=="invalid-label")
    {
      Correlations invalid;
      invalid.addBlock("named").addTrace(Line(forward,"duplicate"),
                                           Line(other,"duplicate"));
      (void)buildTaskGraph(invalid);
      return;
    }
  
  if(argc>1 and std::string(argv[1])=="invalid-block")
    {
      Correlations invalid;
      invalid.blocks.push_back(Correlation{.name="same"});
      invalid.blocks.push_back(Correlation{.name="same"});
      (void)buildTaskGraph(invalid);
      return;
    }
}

void validateScheduling(int argc,char** argv)
{
  const Source backward(0);
  const Source forward(0);
  const CorrelationTaskGraph graph(
    {Task{.id=0,.operation=backward,.inputs{},.labels{}},
     Task{.id=1,.operation=forward,.inputs{},.labels{}},
     Task{.id=2,.operation=ExtendSelectingTPars{Pars(Gamma{.iGamma=1})},
          .inputs{{0,1.0}},.labels{}},
     Task{.id=3,.operation=TraceContraction{},
          .inputs{{2,1.0},{1,1.0}},.labels{}}},
    {OutputBlock{.name="pion",.gammaList{},.traces{3}}});
  
  const CompiledCorrelations stable=
    scheduleCorrelations(graph);
  const CompiledCorrelations outputFirst=
    scheduleCorrelations(graph,scheduleByOutput);
  stable.validate();
  outputFirst.validate();
  
  const auto ids=[](const CompiledCorrelations& compiled)
  {
    std::vector<TaskId> result;
    for(const Task& instruction : compiled.instructions)
      result.push_back(instruction.id);
    
    return result;
  };
  
  assert(ids(stable)==(std::vector<TaskId>{0,1,2,3}));
  assert(ids(outputFirst)==(std::vector<TaskId>{0,2,1,3}));
  assert(outputFirst.outputs==graph.outputs);
  
  if(argc>1 and std::string(argv[1])=="invalid-order")
    {
      (void)scheduleCorrelations(graph,[](const CorrelationTaskGraph&)
      {
        return std::vector<TaskId>{2,0,1,3};
      });
    }
  
  if(argc>1 and std::string(argv[1])=="missing-task")
    {
      (void)scheduleCorrelations(graph,[](const CorrelationTaskGraph&)
      {
        return std::vector<TaskId>{0,2,3};
      });
    }
  
  if(argc>1 and std::string(argv[1])=="invalid-schedule-output")
    {
      CompiledCorrelations invalid=outputFirst;
      invalid.outputs.front().traces.front()=2;
      invalid.validate();
    }
}

void validateLegacyWriter(int argc,char** argv)
{
  const Source source(0);
  const Line backward(Prop{.kappa=0.125,.mass=0.5,.r=1,
                           .charge=0,.mom={1,0,-1},.residue=0.25}*source,"bw");
  const Line forward(Gamma{.iGamma=5}*source,"fw");
  Correlations requests;
  Correlation& block=requests.addBlock("pion");
  block.addGammas(5,5);
  block.addTrace(backward,forward);
  
  const CompiledCorrelations compiled=
    scheduleCorrelations(buildTaskGraph(requests));
  const std::string result=
    writeLegacyInput(compiled);
  
  if(argc>1 and std::string(argv[1])=="dump")
    std::cout<<result;
  
  // Captured from prepareInput's Run printers for the same request
  const std::string expected=
    "NSources 1\n"
    "Name     NoiseType  Tins  Store  \n"
    "SOURCE0  Z4         0     0      \n"
    "\n"
    "NProps 2\n"
    "\n"
    "Name  Ins  SourceName  Tins  Kappa  Mass  R  Charge  ThetaX  ThetaY  ThetaZ  Residue  Store  \n"
    "bw    -    SOURCE0     -1    0.125  0.5   1  0       1       0       -1      0.25     0      \n"
    "fw    G    SOURCE0     -1                 5  0                                        0      \n"
    "\n"
    "NMes2PtsContrBlocks  1\n"
    "\n"
    "NMes2PtsContr 1\n"
    "pion     bw     fw  \n"
    "\n"
    "NGammaContr 1\n"
    "P5P5\n"
    "\n"
    "\n";
  assert(result==expected);
  assert(result==writeLegacyInput(compiled));
  
  const Line weighted(Pars(Gamma{.iGamma=1}),
                      {{2.0,backward},{-0.5,forward}});
  Correlations weightedRequests;
  weightedRequests.addBlock("weighted").addTrace(weighted,forward);
  const std::string weightedOutput=
    writeLegacyInput(scheduleCorrelations(buildTaskGraph(weightedRequests),
                                          scheduleByOutput));
  assert(weightedOutput.find("LINCOMB 2 bw 2 fw -0.5 ")!=std::string::npos);
  
  const Line smeared(Smear{.kappa=0.3,.n=20,.mom={1,2,3}}*source,"smeared");
  const Line phased(Phase{.mom={-1,0,1}}*smeared,"phased");
  Correlations parameterRequests;
  parameterRequests.addBlock("parameters").addTrace(smeared,phased);
  const std::string parameterOutput=
    writeLegacyInput(scheduleCorrelations(buildTaskGraph(parameterRequests)));
  assert(parameterOutput.find("smeared  SM")!=std::string::npos);
  assert(parameterOutput.find("phased   PH")!=std::string::npos);
  assert(parameterOutput.find("0.3")!=std::string::npos);
  assert(parameterOutput.find("20")!=std::string::npos);
  assert(parameterOutput.find("-1")!=std::string::npos);
  
  const Source anotherSource(0);
  const CorrelationTaskGraph nameCollision(
    {Task{.id=0,.operation=source,.inputs{},.labels{}},
     Task{.id=1,.operation=anotherSource,.inputs{},.labels{"SOURCE0"}},
     Task{.id=2,.operation=TraceContraction{},
          .inputs{{0,1.0},{1,1.0}},.labels{}}},
    {OutputBlock{.name="named",.gammaList{},.traces{2}}});
  const std::string namedOutput=
    writeLegacyInput(scheduleCorrelations(nameCollision));
  assert(namedOutput.find("SOURCE1")!=std::string::npos);
  assert(namedOutput.find("named     SOURCE1     SOURCE0")!=std::string::npos);
}

void validateScriptBindings()
{
  Nice nice;
  Correlations& correlations=nice.correlations;
  Evaluator evaluator;
  setupNice(evaluator,correlations);
  registerHostFunction(evaluator,"identity",{"value"},[](const int value)
  {
    return value;
  });
  registerHostFunction(evaluator,"zero",{},[]()
  {
    return 0;
  });
  
  const Source source(0);
  const Line backward(source);
  const Line forward(Pars(Phase{.mom={1,0,0}}),{{1.0,Line(source)}});
  const auto lineInfo=
    evaluator.hostObjectTypes.at(typeid(Line));
  assert(evaluator.hostObjectTypes.at(typeid(Source))->memberAccessorsProvider.contains("tSelect"));
  assert(not evaluator.hostObjectTypes.at(typeid(Source))->memberAccessorsProvider.contains("count"));
  
  const auto value=[](const auto& data)
  {
    return std::make_shared<Value>(Value{.data{data}});
  };
  
  const auto call=[&](const std::string& name,
                      std::vector<std::shared_ptr<Value>> arguments)
  {
    const auto& function=
      std::get<HostFunction>(evaluator.curEnv->at(name).data);
    function(arguments,evaluator);
  };
  
  // Tuple extraction borrows the original host payload, without copying the line
  const std::vector<std::shared_ptr<Value>> arguments{
    value(std::string("pion")),value(HObject{lineInfo,backward}),value(5)};
  const auto concrete=
    hostArgumentsTypeConcretizer<std::string,Line,int>(arguments,"tupleTest");
  assert(&std::get<1>(concrete)==
         std::any_cast<Line>(&std::get<HObject>(arguments[1]->data).data));
  assert(std::apply([](const std::string& name,
                      const Line& line,
                      const int gamma)
  {
    return name=="pion" and not line.isExtender() and gamma==5;
  },concrete));
  assert(std::apply([]()
  {
    return true;
  },hostArgumentsTypeConcretizer<>({},"noArguments")));
  
  call("corrBlock",{value(std::string("pion"))});
  call("corrGammas",{value(std::string("pion")),value(5),value(5)});
  call("corrTrace",{value(std::string("pion")),
                    value(HObject{lineInfo,backward}),
                    value(HObject{lineInfo,forward})});
  
  assert(correlations.blocks.size()==1);
  const Correlation& block=
    correlations.blocks.front();
  assert(block.name=="pion");
  assert(block.gammaList.contains({5,5}));
  assert(block.traces.size()==1);
  assert(block.traces.front()==std::pair(backward,forward));
  
  (*evaluator.curEnv)["bw"]=ValueT{HObject{lineInfo,backward}};
  (*evaluator.curEnv)["fw"]=ValueT{HObject{lineInfo,forward}};
  for(const char* name : {"scriptedSource","scriptedCopy","scriptedOper",
                          "weighted","momentumProduct","legacyText"})
    (*evaluator.curEnv)[name]=std::monostate{};
  
  const auto parser=
    getParser();
  const auto tree=
    createParseTree(parser,R"(
    {
      corrBlock(.name="vector");
      corrGammas(.name="vector",.sink=identity(1),.source=1+zero());
      scriptedSource = Source(0);
      scriptedCopy = scriptedSource;
      sm=Smear(.n=2);
      scriptedOper = Prop(.kappa=0.12,.r=1,.mom=Momentum(1,2,3))
                     * sm * Phase(.mom=Momentum(.z=1))
                     * Gamma(5) * DeltaT(4);
      bw = scriptedOper*scriptedSource;
      fw = Gamma(1)*bw;
      weighted = 2*fw;
      momentumProduct = 2*Momentum(1,2,3);
      corrTrace(.name="vector",.backward=bw,.forward=fw);
      legacyText = writeLegacyInput();
    }
    )");
  
  const CorrelationTaskGraph graph=buildTaskGraph(correlations);
  const CompiledCorrelations compiled=nice.compile(scheduleByOutput);
  assert(compiled.outputs.size()==1);
  assert(compiled.instructions.size()==graph.tasks.size());

  graph.validate();
  const auto executor=
    getParseTreeExecuctor(parser.actionStrings);
  const auto ast=
    executor.execParseTree(tree);
  std::visit(evaluator,*ast);
  
  assert(correlations.blocks.size()==2);
  const Correlation& scripted=
    correlations.blocks.back();
  assert(scripted.name=="vector");
  assert(scripted.gammaList.contains({1,1}));
  const auto& scriptedBackward=
    std::any_cast<const Line&>(std::get<HObject>(evaluator.curEnv->at("bw").data).data);
  const auto& scriptedForward=
    std::any_cast<const Line&>(std::get<HObject>(evaluator.curEnv->at("fw").data).data);
  assert(scripted.traces.front()==std::pair(scriptedBackward,scriptedForward));
  assert(scriptedBackward.isExtender());
  assert(scriptedForward.isExtender());
  const auto& scriptedSource=
    std::any_cast<const Source&>(std::get<HObject>(evaluator.curEnv->at("scriptedSource").data).data);
  const auto& scriptedCopy=
    std::any_cast<const Source&>(std::get<HObject>(evaluator.curEnv->at("scriptedCopy").data).data);
  const auto& scriptedOper=
    std::any_cast<const Oper&>(std::get<HObject>(evaluator.curEnv->at("scriptedOper").data).data);
  const auto& weighted=
    std::any_cast<const Line&>(std::get<HObject>(evaluator.curEnv->at("weighted").data).data);
  const auto& momentumProduct=
    std::any_cast<const Momentum&>(std::get<HObject>(evaluator.curEnv->at("momentumProduct").data).data);
  assert(scriptedSource==scriptedCopy);
  assert(scriptedBackward==scriptedOper*scriptedSource);
  assert(weighted==2*scriptedForward);
  assert(momentumProduct==Momentum({2,4,6}));
  const std::string& legacyText=
    std::get<std::string>(evaluator.curEnv->at("legacyText").data);
  assert(legacyText==writeLegacyInput(nice.compile()));
  assert(legacyText.find("NMes2PtsContrBlocks  2\n")!=std::string::npos);
  assert(legacyText.find("NGammaContr 1\nV1V1\n")!=std::string::npos);
}

/// Small parameter type demonstrating the registration contract
struct ExamplePars
{
  int t{4};
  double weight{1.0};
  std::string label{"source"};
  bool enabled{true};

  PROVIDE_NICE_MEMBERS(ExamplePars,
		 MEMB(t),
		 MEMB(weight),
		 MEMB(label),
		 MEMB(enabled))
};

/// Concrete nested types; HObject remains an interpreter implementation detail
struct ExampleGroup
{
  ExamplePars pars;
  int tag{1};
};

struct ExampleEnvelope
{
  ExampleGroup group;
};

/// Runs through the real parser and evaluator, then checks the concrete C++ values
void validateHostObjects(const std::string& script)
{
  Evaluator ev;
  const auto parsInfo=REGISTER_HOBJECT(ev,ExamplePars);
  const auto groupInfo=registerHObject<ExampleGroup>(ev,"ExampleGroup",
    std::make_tuple(HostMemberInfo{"pars",&ExampleGroup::pars,parsInfo},
                    HostMemberInfo{"tag",&ExampleGroup::tag}));
  registerHObject<ExampleEnvelope>(ev,"ExampleEnvelope",
    std::make_tuple(HostMemberInfo{"group",&ExampleEnvelope::group,groupInfo}));
  // Keep the checked variables in the outer scope of the script block
  for(const char* name : {"p","q","defaults","temporary","nested","copy","detached","nestedTemporary","replaced"})
    (*ev.curEnv)[name]=std::monostate{};
  const auto parser=getParser();
  const auto tree=createParseTree(parser,script.c_str());
  const auto executor=getParseTreeExecuctor(parser.actionStrings);
  const auto ast=executor.execParseTree(tree);
  std::visit(ev,*ast);

  auto pars=[&](const std::string& name)->const ExamplePars&
  {
    return std::any_cast<const ExamplePars&>(std::get<HObject>(ev.curEnv->at(name).data).data);
  };
  const auto& p=pars("p");
  const auto& q=pars("q");
  const auto& defaults=pars("defaults");
  if(p.t!=7 or p.weight!=2.5 or p.label!="original" or not p.enabled or
     q.t!=10 or q.weight!=5.0 or q.label!="copy" or q.enabled or
     defaults.t!=4 or defaults.weight!=1.0 or defaults.label!="source" or not defaults.enabled or
     std::get<int>(ev.curEnv->at("temporary").data)!=9)
    errorEmitter("host object registration test failed");
  auto envelope=[&](const std::string& name)->const ExampleEnvelope&
  {
    return std::any_cast<const ExampleEnvelope&>(std::get<HObject>(ev.curEnv->at(name).data).data);
  };
  const auto& nested=envelope("nested");
  const auto& copy=envelope("copy");
  if(nested.group.pars.t!=14 or nested.group.pars.weight!=6.0 or nested.group.tag!=17 or
     copy.group.pars.t!=7 or copy.group.pars.weight!=2.5 or
     pars("detached").t!=99 or envelope("replaced").group.pars.t!=21 or
     envelope("replaced").group.pars.weight!=1.0 or
     std::get<int>(ev.curEnv->at("nestedTemporary").data)!=4)
    errorEmitter("nested host reference test failed");
}

int main(int argc,char** argv)
{
  validateCorrParts();
  validateTaskGraph(argc,argv);
  validateGraphCompiler(argc,argv);
  validateScheduling(argc,argv);
  validateLegacyWriter(argc,argv);
  validateScriptBindings();
  validateHostObjects(R"(
    {
      defaults = ExamplePars();
      p = ExamplePars(7, .label="original", .weight=2.5);
      q = p;
      q.t += 3;
      q.weight *= 2;
      q.label = "copy";
      q.enabled = 0;
      temporary = ExamplePars(.t=9).t;
      nested = ExampleEnvelope(.group=ExampleGroup(.pars=p));
      copy = nested;
      detached = nested.group.pars;
      detached.t = 99;
      nested.group.pars.t += 5;
      nested.group.pars.weight = 6;
      nested.group.tag = 8;
      fun changeSibling()
      {
        nested.group.tag = 17;
        return 2;
      }
      nested.group.pars.t += changeSibling();
      replaced = nested;
      replaced.group.pars = ExamplePars(.t=21);
      nestedTemporary = ExampleEnvelope().group.pars.t;
    }
    )");
}
