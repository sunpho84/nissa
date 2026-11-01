#include "corrParts.hpp"
#include "correlations.hpp"
#include "correlationTaskGraph.hpp"
#include "graphCompiler.hpp"
#include "schedulingStrategies.hpp"
#include "legacyInputWriter.hpp"
#include "niceSetup.hpp"
#include "nice.hpp"
#include "nASTy.hpp"
#include "new_types/dirac.hpp"
#include "routines/mpi_routines.hpp"

#include <algorithm>
#include <cassert>
#include <complex>
#include <cstdio>
#include <unistd.h>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace nissa::nice;

/// Find a gamma extension carrying the requested time selection in a line tree
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

/// Check aggregate defaults, source identities, adjoints, and weighted line construction
void validateCorrParts()
{
  // Propagators inherit every default while explicit members override selected values
  static_assert(std::is_aggregate_v<Prop>);
  const Prop initialDefaults=Prop::defaults;
  assert(initialDefaults.kappa==0.125 and initialDefaults.r==0);
  Prop::defaults=Prop{.kappa=0.13,.mass=0.25,.r=-1,.charge=0.5,
                      .mom={1,2,3},.residue=1e-12};
  const Prop inherited;
  assert(inherited==Prop::defaults);
  const Prop overridden{.kappa=0.14,.r=1};
  assert(overridden.kappa==0.14 and overridden.r==1);
  assert(overridden.mass==0.25 and overridden.mom==SpatMomentum({1,2,3}));
  Prop::defaults=initialDefaults;
  assert(inherited.r==-1 and inherited.kappa==0.13);
  
  // Smearing defaults obey the same snapshot rule and preserve aggregate initialization
  static_assert(std::is_aggregate_v<Smear>);
  const Smear initialSmearDefaults=Smear::defaults;
  assert(initialSmearDefaults==Smear({.kappa=0,.n=0,.mom={}}));
  Smear::defaults=Smear{.kappa=0.3,.n=20,.mom={1,2,3}};
  const Smear inheritedSmear;
  assert(inheritedSmear==Smear::defaults);
  assert(Smear{.n=0}.kappa==0.3 and Smear{.n=0}.mom==SpatMomentum({1,2,3}));
  Smear::defaults=initialSmearDefaults;
  assert(inheritedSmear.n==20);
  
  // Copies retain source identity, but a fresh source at the same time is independent
  const Source eta(0);
  const Source copy=
    eta;
  const Source independent(0);
  assert(copy==eta);
  assert(independent!=eta);
  
  // Propagator conjugation reverses the Wilson sign and is an involution
  const Prop prop{.kappa=0.1394267,.mass=0.00066690,.r=-1,
                  .residue=1e-20};
  const Prop reversed=
    prop.dag();
  assert(reversed.r==1);
  assert(reversed.dag()==prop);
  
  // Momentum scaling and phase conjugation preserve the spatial-component convention
  const SpatMomentum momentum{1,2,3};
  assert(momentum[0]==1 and momentum[2]==3);
  assert(((2.0*momentum).asArray()==std::array<double,3>{2,4,6}));
  assert(Phase{.mom=momentum}.dag().mom==SpatMomentum({-1,-2,-3}));
  
  // Product adjoints and scalar weights survive copying and a second adjoint
  const Oper op=
    prop*Phase{.mom=momentum}*Gamma{.iGamma=5};
  const Oper reversedOp=
    op.dag();
  assert(reversedOp.dag()==op);
  assert((2.0*op).weight==2.0);
  
  // Display labels do not change line identity or erase selected times
  const Line original=
    op*DeltaT{.t=5}*eta;
  const Line named(original,"first");
  const Line renamed(original,"second");
  assert(named==renamed);
  assert(named.describe()==renamed.describe());
  assert(containsSelection(original,5));
  
  // An explicitly unmerged sum retains separate weighted input terms
  const Line sum=
    original+2.0*original;
  const auto* extender=
    std::get_if<Extender>(&sum.data);
  assert(extender and extender->rhs.size()==2);
  assert(extender->rhs[0].first==1.0);
  assert(extender->rhs[1].first==2.0);
  assert(addLines(original,original,false)!=original);
}

/// Check adjoints against NISSA's actual matrices, including product reversal and signs
void validateGammaAdjoints()
{
  // Use the backend matrices as the reference rather than duplicating an expected-sign table
  nissa::init_base_gamma();
  using Matrix=
    std::array<std::array<std::complex<double>,4>,4>;
  std::array<Matrix,16> matrices{};
  
  // Expand the backend sparse representation: each row stores one nonzero entry
  for(std::size_t index=0;index<matrices.size();index++)
    for(std::size_t row=0;row<4;row++)
      {
        const auto& gamma=nissa::base_gamma[index];
        matrices[index][row][gamma.pos[row]]={gamma.entr[row][0],gamma.entr[row][1]};
      }
  
  // Dense multiplication provides an independent check of symbolic product order
  const auto multiply=[](const Matrix& left,
                         const Matrix& right)
  {
    Matrix result{};
    for(std::size_t row=0;row<4;row++)
      for(std::size_t column=0;column<4;column++)
        for(std::size_t inner=0;inner<4;inner++)
          result[row][column]+=left[row][inner]*right[inner][column];
    
    return result;
  };
  
  // The reference adjoint is the complex-conjugate transpose
  const auto adjoint=[](const Matrix& matrix)
  {
    Matrix result{};
    for(std::size_t row=0;row<4;row++)
      for(std::size_t column=0;column<4;column++)
        result[row][column]=std::conj(matrix[column][row]);
    
    return result;
  };
  
  // Evaluate a symbolic gamma product numerically, including every node weight
  const auto matrixOf=[&](auto&& self,
                         const Oper& operation)->Matrix
  {
    Matrix result=matrices[std::get<Gamma>(operation.pars).iGamma];
    
    if(operation.rhs)
      result=multiply(result,self(self,*operation.rhs));
    
    for(auto& row : result)
      for(auto& entry : row)
        entry*=operation.weight;
    
    return result;
  };
  
  // Entries and weights are exactly representable here, so exact matrix equality is sufficient
  for(std::size_t left=0;left<matrices.size();left++)
    {
      // Check both the direct gamma API and dispatch through the parameter variant
      const Gamma gamma{.iGamma=left};
      assert(matrixOf(matrixOf,gamma.dag())==adjoint(matrices[left]));
      assert(dag(Pars{gamma})==gamma.dag());
      
      // All ordered pairs test product reversal, accumulated signs, and a nonunit real weight
      for(std::size_t right=0;right<matrices.size();right++)
        {
          const Oper product=2.0*(Oper(gamma)*Oper(Gamma{.iGamma=right}));
          assert(matrixOf(matrixOf,product.dag())==adjoint(matrixOf(matrixOf,product)));
        }
    }
  
  const Source source(0);
  
  // The anti-Hermitian sign must also reach the coefficient of a constructed line
  const Line signedLine=Gamma{.iGamma=6}.dag()*source;
  const auto& extension=std::get<Extender>(signedLine.data);
  assert(extension.rhs.front().first==-1.0);
}

/// Check graph invariants and parameter ranges, with named modes for rejected inputs
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
  
  // All standard gamma entries remain valid through scheduling and serialization
  for(std::size_t index=0;index<16;index++)
    {
      auto& extension=
        std::get<ExtendSelectingTPars>(graph.tasks[1].operation);
      extension.pars=Gamma{.iGamma=index};
      graph.validate();
      const auto compiled=
        scheduleCorrelations(graph);
      assert(compiled.instructions.size()==2);
      assert(writeLegacyInput(compiled).find("NProps 1")!=std::string::npos);
    }
  
  std::get<ExtendSelectingTPars>(graph.tasks[1].operation).pars=Gamma{.iGamma=1};
  
  // Reject the first unsupported index, historical extra entries, and wrapped negative values
  const std::pair<const char*,std::size_t> invalidGammaChecks[]{
    {"invalid-gamma-boundary",16},
    {"invalid-gamma-extra-17",17},
    {"invalid-gamma-extra-18",18},
    {"invalid-gamma-large",999},
    {"invalid-gamma-wrapped",std::numeric_limits<std::size_t>::max()}};
  
  for(const auto& [mode,index] : invalidGammaChecks)
    if(argc>1 and std::string(argv[1])==mode)
      {
        std::get<ExtendSelectingTPars>(graph.tasks[1].operation).pars=Gamma{.iGamma=index};
        graph.validate();
        errorEmitter("invalid extension gamma was unexpectedly accepted");
      }
  
  // Validation must also catch an invalid parameter introduced after scheduling
  if(argc>1 and std::string(argv[1])=="invalid-compiled-gamma")
    {
      auto compiled=
        scheduleCorrelations(graph);
      std::get<ExtendSelectingTPars>(compiled.instructions[1].operation).pars=Gamma{.iGamma=999};
      writeLegacyInput(compiled);
      errorEmitter("invalid compiled gamma was unexpectedly serialized");
    }
  
  // The initial Wilson default is unset; both physical signs remain valid
  assert(Prop{}.r==0);
  
  for(const int sign : {-1,1})
    {
      const Prop prop{.r=sign};
      assert(prop.dag().r==-sign);
      std::get<ExtendSelectingTPars>(graph.tasks[1].operation).pars=prop;
      graph.validate();
      assert(not writeLegacyInput(scheduleCorrelations(graph)).empty());
    }
  
  // Unset and nonphysical signs must fail before the legacy sign mapping
  const std::pair<const char*,int> invalidWilsonChecks[]{
    {"unset-wilson-sign",0},
    {"invalid-wilson-positive",7},
    {"invalid-wilson-negative",-7}};
  
  for(const auto& [mode,sign] : invalidWilsonChecks)
    if(argc>1 and std::string(argv[1])==mode)
      {
        std::get<ExtendSelectingTPars>(graph.tasks[1].operation).pars=Prop{.r=sign};
        graph.validate();
        errorEmitter("invalid Wilson sign was unexpectedly accepted");
      }
  
  if(argc>1 and std::string(argv[1])=="invalid-compiled-wilson-sign")
    {
      auto compiled=
        scheduleCorrelations(graph);
      std::get<ExtendSelectingTPars>(compiled.instructions[1].operation).pars=Prop{.r=7};
      writeLegacyInput(compiled);
      errorEmitter("invalid compiled Wilson sign was unexpectedly serialized");
    }
  
  if(argc>1 and std::string(argv[1])=="unset-wilson-dagger")
    {
      Prop{}.dag();
      errorEmitter("unset Wilson sign was unexpectedly conjugated");
    }
  
  // Accept the iteration-count endpoints without performing any numerical smearing
  for(const double count : {0.0,static_cast<double>(std::numeric_limits<int>::max())})
    {
      auto smearingGraph=graph;
      std::get<ExtendSelectingTPars>(smearingGraph.tasks[1].operation).pars=Smear{.n=count};
      smearingGraph.validate();
      assert(not writeLegacyInput(scheduleCorrelations(smearingGraph)).empty());
    }
  
  // Fractional, negative, nonfinite, and excessive counts must not reach integer serialization
  const std::pair<const char*,double> invalidSmearChecks[]{
    {"fractional-smear-count",2.5},
    {"negative-smear-count",-1},
    {"overflow-smear-count",static_cast<double>(std::numeric_limits<int>::max())+1},
    {"nan-smear-count",std::numeric_limits<double>::quiet_NaN()},
    {"infinite-smear-count",std::numeric_limits<double>::infinity()}};
  
  for(const auto& [mode,count] : invalidSmearChecks)
    if(argc>1 and std::string(argv[1])==mode)
      {
        std::get<ExtendSelectingTPars>(graph.tasks[1].operation).pars=Smear{.n=count};
        graph.validate();
        errorEmitter("invalid smearing count was unexpectedly accepted");
      }
  
  if(argc>1 and std::string(argv[1])=="invalid-compiled-smear-count")
    {
      auto compiled=scheduleCorrelations(graph);
      std::get<ExtendSelectingTPars>(compiled.instructions[1].operation).pars=Smear{.n=2.5};
      writeLegacyInput(compiled);
      errorEmitter("invalid compiled smearing count was unexpectedly serialized");
    }
  
  // Sources and line extensions produce fields that may be contracted
  const CorrelationTaskGraph fieldInputs(
    {Task{.id=0,.operation=source,.inputs{},.labels{}},
     Task{.id=1,.operation=ExtendSelectingTPars{Pars(Gamma{.iGamma=1})},
          .inputs{{0,1.0}},.labels{}},
     Task{.id=2,.operation=TraceContraction{},.inputs{{0,1.0},{1,1.0}},.labels{}}},
    {OutputBlock{.name="fields",.gammaList{},.traces{2}}});
  const auto validFields=scheduleCorrelations(fieldInputs);
  assert(validFields.instructions.size()==3);
  assert(not writeLegacyInput(validFields).empty());
  
  // Contraction results cannot feed an extension or either operand of another trace
  const std::pair<const char*,Task> invalidFieldInputs[]{
    {"extension-from-trace",Task{.id=3,
      .operation=ExtendSelectingTPars{Pars(Gamma{.iGamma=1})},
      .inputs{{2,1.0}},.labels{}}},
    {"weighted-extension-from-trace",Task{.id=3,
      .operation=ExtendSelectingTPars{Pars(Gamma{.iGamma=1})},
      .inputs{{0,1.0},{2,0.5}},.labels{}}},
    {"trace-from-trace-backward",Task{.id=3,.operation=TraceContraction{},
      .inputs{{2,1.0},{0,1.0}},.labels{}}},
    {"trace-from-trace-forward",Task{.id=3,.operation=TraceContraction{},
      .inputs{{0,1.0},{2,1.0}},.labels{}}}};
  
  for(const auto& [mode,task] : invalidFieldInputs)
    if(argc>1 and std::string(argv[1])==mode)
      {
        auto invalid=fieldInputs;
        invalid.tasks.push_back(task);
        invalid.validate();
        errorEmitter("contraction was unexpectedly accepted as a field input");
      }
  
  if(argc>1 and std::string(argv[1])=="compiled-extension-from-trace")
    {
      auto invalid=validFields;
      invalid.instructions.push_back(invalidFieldInputs[0].second);
      invalid.validate();
      errorEmitter("compiled contraction was unexpectedly accepted as a field input");
    }
  
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

/// Check stable task identities, shared dependencies, coefficients, and output ordering
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

/// Check complete dependency-compatible schedules, including unused and independent tasks
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
  
  // A schedule must retain unreferenced tasks and may reorder independent ones
  const CorrelationTaskGraph independent(
    {Task{.id=0,.operation=backward,.inputs{},.labels{}},
     Task{.id=1,.operation=forward,.inputs{},.labels{}}},{});
  assert(ids(scheduleCorrelations(independent))==(std::vector<TaskId>{0,1}));
  assert(ids(scheduleCorrelations(independent,scheduleByOutput))==(std::vector<TaskId>{0,1}));
  assert(ids(scheduleCorrelations(independent,[](const CorrelationTaskGraph&)
  {
    return std::vector<TaskId>{1,0};
  }))==(std::vector<TaskId>{1,0}));
  
  // Both built-in strategies accept an empty graph
  const CorrelationTaskGraph empty({},{});
  assert(scheduleCorrelations(empty).instructions.empty());
  assert(scheduleCorrelations(empty,scheduleByOutput).instructions.empty());
  
  // A complete schedule needs the right count, unique IDs, and valid references
  const std::pair<const char*,std::vector<TaskId>> invalidSchedules[]{
    {"omitted-trailing-task",{0}},
    {"excess-scheduled-tasks",{0,1,0}},
    {"duplicate-scheduled-task",{0,0}},
    {"nonexistent-scheduled-task",{0,2}}};
  
  for(const auto& [mode,order] : invalidSchedules)
    if(argc>1 and std::string(argv[1])==mode)
      {
        scheduleCorrelations(independent,[&order](const CorrelationTaskGraph&)
        {
          return order;
        });
        errorEmitter("invalid schedule was unexpectedly accepted");
      }
  
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

/// Compare legacy tables with a reference fixture and reject names interpreted as syntax
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
  // Near-matches remain ordinary legacy tokens unless the reader treats them as syntax
  for(const char* name : {"LINCOMB_1","xLINCOMB","x/*comment","*/tail","//name"})
    {
      auto safe=compiled;
      const auto found=
        std::find_if(safe.instructions.begin(),safe.instructions.end(),[](const Task& task)
        {
          return std::holds_alternative<Source>(task.operation);
        });
      assert(found!=safe.instructions.end());
      found->labels={name};
      assert(writeLegacyInput(safe).find(name)!=std::string::npos);
    }
  
  // Exercise source labels (0), line labels (1), and output block names (2)
  const std::tuple<const char*,const char*,int> invalidNames[]{
    {"reserved-source-name","LINCOMB",0},
    {"reserved-line-name","lincomb",1},
    {"reserved-block-name","LiNcOmB",2},
    {"comment-source-name","/*hidden",0},
    {"comment-line-name","/*",1},
    {"comment-block-name","/**/",2},
    {"whitespace-source-name","bad name",0}};
  
  for(const auto& [mode,name,role] : invalidNames)
    if(argc>1 and std::string(argv[1])==mode)
      {
        auto invalid=compiled;
        
        if(role==2)
          invalid.outputs.front().name=name;
        else
          {
            const auto found=
              std::find_if(invalid.instructions.begin(),invalid.instructions.end(),[role](const Task& task)
              {
                if(role==0)
                  return std::holds_alternative<Source>(task.operation);
                else
                  return std::holds_alternative<ExtendSelectingTPars>(task.operation);
              });
            assert(found!=invalid.instructions.end());
            found->labels={name};
          }
        
        writeLegacyInput(invalid);
        errorEmitter("legacy syntax was unexpectedly accepted as a name");
      }

}

/// Exercise registered NICE functions and construction-time defaults through the script evaluator
void validateScriptBindings(int argc,
                            char** argv)
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
  
  if(argc>1 and std::string(argv[1])=="invalid-default-wilson-sign")
    {
      call("setDefaultPropR",{value(7)});
      errorEmitter("invalid default Wilson sign was unexpectedly accepted");
    }
  
  if(argc>1 and std::string(argv[1])=="invalid-default-smear-count")
    {
      call("setDefaultSmearN",{value(-1)});
      errorEmitter("invalid default smearing count was unexpectedly accepted");
    }
  
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
                          "weighted","sum","momentumProduct","legacyText"})
    (*evaluator.curEnv)[name]=std::monostate{};
  
  const auto parser=
    getParser();
  // Capture constructed values to distinguish inheritance, explicit overrides, and later default changes
  std::vector<Smear> observedSmears;
  registerHostFunction(evaluator,"recordSmear",{"smear"},[&observedSmears](const Smear& smear)
  {
    observedSmears.push_back(smear);
  });
  std::vector<Prop> observedProps;
  registerHostFunction(evaluator,"recordProp",{"prop"},[&observedProps](const Prop& prop)
  {
    observedProps.push_back(prop);
  });
  const auto defaultsTree=
    createParseTree(parser,R"(
      setDefaultPropR(-1);
      first=Prop();
      recordProp(first);
      recordProp(Prop(.kappa=0.13,.r=1));
      setDefaultPropR(1);
      recordProp(first);
      recordProp(Prop());
      recordProp(Prop(.r=0));
      setDefaultPropR(0);
      recordProp(Prop());
      setDefaultPropKappa(0.14);
      setDefaultPropMass(0.02);
      recordProp(Prop());
      recordProp(Prop(.kappa=0.13,.mass=0));
      recordProp(first);
      setDefaultPropKappa(0.125);
      setDefaultPropMass(0.0);
      recordSmear(Smear());
      setDefaultSmearKappa(0.3);
      setDefaultSmearN(20);
      smearing=Smear();
      recordSmear(smearing);
      recordSmear(Smear(.kappa=0.4,.n=0));
      setDefaultSmearKappa(0.5);
      setDefaultSmearN(30);
      recordSmear(smearing);
      recordSmear(Smear());
      setDefaultSmearKappa(0.0);
      setDefaultSmearN(0);
      recordSmear(Smear());
    )");
  const auto defaultsExecutor=
    getParseTreeExecuctor(parser.actionStrings);
  const auto defaultsAst=
    defaultsExecutor.execParseTree(defaultsTree);
  std::visit(evaluator,*defaultsAst);
  assert(observedProps.size()==9);
  assert(observedProps[0].r==-1 and observedProps[0].kappa==0.125);
  assert(observedProps[1].r==1 and observedProps[1].kappa==0.13);
  assert(observedProps[2]==observedProps[0]);
  assert(observedProps[3].r==1 and observedProps[3].kappa==0.125);
  assert(observedProps[4].r==0 and observedProps[5].r==0);
  assert(observedProps[6].kappa==0.14 and observedProps[6].mass==0.02);
  assert(observedProps[7].kappa==0.13 and observedProps[7].mass==0);
  assert(observedProps[8]==observedProps[0]);
  assert(observedSmears.size()==6);
  assert(observedSmears[0].kappa==0 and observedSmears[0].n==0);
  assert(observedSmears[1].kappa==0.3 and observedSmears[1].n==20);
  assert(observedSmears[2].kappa==0.4 and observedSmears[2].n==0);
  assert(observedSmears[3]==observedSmears[1]);
  assert(observedSmears[4].kappa==0.5 and observedSmears[4].n==30);
  assert(observedSmears[5]==observedSmears[0]);
  assert(Prop{}.r==0);
  
  const auto tree=
    createParseTree(parser,R"(
    {
      corrBlock(.name="vector");
      corrGammas(.name="vector",.sink=identity(1),.source=1+zero());
      scriptedSource = Source(0);
      scriptedCopy = scriptedSource;
      sm=Smear(.n=2);
      scriptedOper = Prop(.kappa=0.12,.r=1,.mom=SpatMomentum(1,2,3))
                     * sm * Phase(.mom=SpatMomentum(.z=1))
                     * Gamma(5) * DeltaT(4);
      bw = scriptedOper*scriptedSource;
      fw = Gamma(1)*bw;
      weighted = fw;
      weighted *= 2;
      sum = bw;
      sum += weighted;
      momentumProduct = 2*SpatMomentum(1,2,3);
      momentumProduct /= 2;
      momentumProduct = momentumProduct/0.5;
      corrTrace(.name="vector",.backward=bw,.forward=fw);
      corrTrace(.name="vector",.backward=sum,.forward=fw);
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
  const auto& sum=
    std::any_cast<const Line&>(std::get<HObject>(evaluator.curEnv->at("sum").data).data);
  const auto& momentumProduct=
    std::any_cast<const SpatMomentum&>(std::get<HObject>(evaluator.curEnv->at("momentumProduct").data).data);
  assert(scriptedSource==scriptedCopy);
  assert(scriptedBackward==scriptedOper*scriptedSource);
  assert(weighted==2*scriptedForward);
  assert(sum==scriptedBackward+weighted);
  assert(scripted.traces.size()==2);
  assert(scripted.traces.back()==std::pair(sum,scriptedForward));
  assert(momentumProduct==SpatMomentum({2,4,6}));
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

/// Host type with binary subtraction and division but no compound operators
struct ExampleNumber
{
  /// Integer payload exposed as a script member
  int value{};
  
  /// Return a new value after subtraction
  ExampleNumber operator-(const int rhs) const
  {
    return {value-rhs};
  }
  
  /// Return a new value after division
  ExampleNumber operator/(const int rhs) const
  {
    return {value/rhs};
  }
};

/// Check program scope, expression evaluation, returns, and reference arguments
void validateProgramGrammar(int argc,
                            char** argv)
{
  const auto parser=getParser();
  const auto executor=getParseTreeExecuctor(parser.actionStrings);
  
  int observed{};
  Evaluator evaluator;
  registerHostFunction(evaluator,"record",{"value"},[&observed](const int value)
  {
    observed=value;
  });
  std::string observedString;
  registerHostFunction(evaluator,"recordString",{"value"},[&observedString](const std::string& value)
  {
    observedString=value;
  });
  registerHostFunction(evaluator,"minimumInteger",{},[]()
  {
    return std::numeric_limits<int>::min();
  });
  const auto numberInfo=
    registerHObject<ExampleNumber>(evaluator,"ExampleNumber",
      std::make_tuple(HostMemberInfo{"value",&ExampleNumber::value}));
  registerHObjectBinaryOperator<DifferenceOperator,ExampleNumber,int>(numberInfo,nullptr,numberInfo);
  registerHObjectBinaryOperator<DivisionOperator,ExampleNumber,int>(numberInfo,nullptr,numberInfo);
  
  const auto unbracedTree=createParseTree(parser,"value=1; record(value+2);");
  const auto unbraced=executor.execParseTree(unbracedTree);
  const auto* unbracedBody=std::get_if<ASTNodesNode>(unbraced.get());
  assert(unbracedBody and unbracedBody->list.size()==2);
  std::visit(evaluator,*unbraced);
  assert(observed==3);
  
  const auto bracedTree=createParseTree(parser,"{ value=4; record(value+2); }");
  const auto braced=executor.execParseTree(bracedTree);
  const auto* bracedBody=std::get_if<ASTNodesNode>(braced.get());
  assert(bracedBody and bracedBody->list.size()==1);
  const auto* block=std::get_if<ASTNodesNode>(bracedBody->list.front().get());
  assert(block and block->list.size()==2);
  std::visit(evaluator,*braced);
  assert(observed==6);
  
  const auto emptyTree=createParseTree(parser,"");
  const auto empty=executor.execParseTree(emptyTree);
  const auto* emptyBody=std::get_if<ASTNodesNode>(empty.get());
  assert(emptyBody and emptyBody->list.empty());
  std::visit(evaluator,*empty);
  
  // String tokens retain escaped quotes and backslashes until decoding
  const std::pair<const char*,std::string> stringChecks[]{
    {R"(recordString("");)",""},
    {R"(recordString("plain");)","plain"},
    {R"(recordString("a\"b");)","a\"b"},
    {R"(recordString("\\");)","\\"},
    {R"(recordString("\\\"");)","\\\""},
    {R"(recordString("left\"middle\"right");)","left\"middle\"right"},
    {R"(recordString("a\nb\tc\rd\be\f");)","a\nb\tc\rd\be\f"},
    {R"(recordString("// not a comment /* either */");)","// not a comment /* either */"},
    {R"(recordString("a\""+"b");)","a\"b"},
    {R"(recordString("one");recordString("two\"");)","two\""}};
  
  for(const auto& [script,expected] : stringChecks)
    {
      observedString="not recorded";
      const auto tree=
        createParseTree(parser,script);
      const auto ast=
        executor.execParseTree(tree);
      std::visit(evaluator,*ast);
      assert(observedString==expected);
    }
  
  if(argc>1 and std::string(argv[1])=="unterminated-script-string")
    {
      createParseTree(parser,R"(recordString("a\");)");
      errorEmitter("unterminated script string was unexpectedly accepted");
    }
  
  // Comments must compose with tokens, division, CRLF, and end-of-input
  const char* commentChecks[]{
    "/**/ record(21);",
    "/***/record(21);",
    "/****/record(21);",
    "/* stars **inside and / slash */record(21);",
    "/* first\nsecond */record(21);",
    "//\nrecord(21);",
    "//\r\nrecord(21);",
    "// comment\rrecord(21);",
    "record(21);// trailing",
    "record(21);//",
    "\r\nrecord(\r\n21\r\n);\r\n",
    "record(84/**//4);",
    "record(20/* first */+1/* second */);"};
  
  for(const char* script : commentChecks)
    {
      observed=0;
      const auto tree=
        createParseTree(parser,script);
      const auto ast=
        executor.execParseTree(tree);
      std::visit(evaluator,*ast);
      assert(observed==21);
    }
  
  for(const char* script : {"//","// trailing","/**/","/***/"," \r\n\t"})
    {
      const auto tree=
        createParseTree(parser,script);
      const auto ast=
        executor.execParseTree(tree);
      const auto* body=
        std::get_if<ASTNodesNode>(ast.get());
      assert(body and body->list.empty());
    }
  
  if(argc>1 and std::string(argv[1])=="unterminated-script-comment")
    {
      createParseTree(parser,"record(21); /* missing end");
      errorEmitter("unterminated script comment was unexpectedly accepted");
    }
  
  // Mixed operators follow the shared C/C++ precedence and associativity groups
  const std::pair<const char*,int> precedenceChecks[]{
    {"record(10-3-2);",5},
    {"record(16/4/2);",2},
    {"record(8-3+2);",7},
    {"record(8*3/2);",12},
    {"record(8%3*2);",4},
    {"record(20/3%2);",0},
    {"record(2+3*4);",14},
    {"record((2+3)*4);",20},
    {"record((3>2<1)+0);",0},
    {"record((0==1!=2)+0);",1},
    {"record((0==1<2)+0);",0},
    {"record((1 or 0 and 0)+0);",1},
    {"a=b=7; record(a*10+b);",77},
    {"a=1; b=2; c=3; a=b=c; record(a*100+b*10+c);",333}};
  
  for(const auto& [script,expected] : precedenceChecks)
    {
      const auto tree=
        createParseTree(parser,script);
      const auto ast=
        executor.execParseTree(tree);
      std::visit(evaluator,*ast);
      assert(observed==expected);
    }
  
  // Short-circuited operands must not execute side effects or invalid expressions
  const std::pair<const char*,int> logicalChecks[]{
    {"x=0; 0 && x++; record(x);",0},
    {"x=0; 1 || x++; record(x);",0},
    {"x=0; 0 and x++; record(x);",0},
    {"x=0; 1 or x++; record(x);",0},
    {"x=0; y=1 && x++; record(y*10+x);",1},
    {"x=0; y=0 || x++; record(y*10+x);",1},
    {"x=0; fun step(){x=x+1; return 2;} y=1 and step(); record(y*10+x);",11},
    {"x=0; fun step(){x=x+1; return 2;} y=0 or step(); record(y*10+x);",11},
    {"record((0 && (1/0))+0);",0},
    {"record((1 || missing())+0);",1},
    {"record(((0 and missing()) or 2)+0);",1},
    {"x=0; y=(x++==0) && (x==1); record(y*10+x);",11},
    {"record((0.5 and -2)+0);",1},
    {"record((0.0 or 0)+0);",0}};
  
  for(const auto& [script,expected] : logicalChecks)
    {
      const auto tree=
        createParseTree(parser,script);
      const auto ast=
        executor.execParseTree(tree);
      observed=-1;
      std::visit(evaluator,*ast);
      assert(observed==expected);
    }
  
  // Returns propagate through nested control flow and stop at function-call boundaries
  const std::pair<const char*,int> returnChecks[]{
    {"fun f(){if(1) return 7; return 9;} record(f());",7},
    {"fun f(){if(0) return 7; else return 8; return 9;} record(f());",8},
    {"fun f(){{return 6;} return 9;} record(f());",6},
    {"fun f(){for(i=0;i<3;i++){if(i==1) return 4;} return 9;} record(f());",4},
    {"i=0; fun f(){for(i=0;i<3;i++){return 5;} return 9;} x=f(); record(x*10+i);",50},
    {"fun f(){for(i=0;i<2;i++){for(j=0;j<2;j++){return 6;}} return 9;} record(f());",6},
    {"fun inner(){return 3;} fun outer(){inner(); return 7;} record(outer());",7},
    {"fun f(){if(1) 42; 12; return 6;} record(f());",6},
    {"fun f(){a=13; return a;} record(f());",13},
    {"fun f(){return 0;} record(f());",0},
    {"fun f(){a=1;} f(); record(6);",6}};
  
  for(const auto& [script,expected] : returnChecks)
    {
      const auto tree=
        createParseTree(parser,script);
      const auto ast=
        executor.execParseTree(tree);
      observed=-1;
      std::visit(evaluator,*ast);
      assert(observed==expected);
    }
  
  // Compound assignment preserves native conversions and updates the original value cell
  const std::pair<const char*,int> compoundChecks[]{
    {"x=1; x+=2; record(x);",3},
    {"x=9; x-=4; record(x);",5},
    {"x=3; x*=4; record(x);",12},
    {"x=9; x/=2; record(x);",4},
    {"x=3; x+=0.75; record(x);",3},
    {"x=5; x/=0.5; record(x);",10},
    {"x=1.0; x+=2; record((x==3.0)+0);",1},
    {"x=2; x+=x; record(x);",4},
    {"x=2; record((x*=3)+x);",12},
    {"a=makeVec(2,3); a[0]+=a[1]; record(a[0]);",5},
    {"struct S{x;}; s=S(4); s.x*=3; record(s.x);",12},
    {"fun inc(&a){a+=2;} x=1; inc(&x); record(x);",3},
    {"x=1; fun rhs(){x=10; return 2;} x+=rhs(); record(x);",12},
    {"i=0; a=makeVec(1,2); a[i++]+=3; record(i*10+a[0]);",14},
    {"x=\"a\"; x+=\"b\"; record((x==\"ab\")+0);",1},
    {"x=1; x=\"hello\"; record((x==\"hello\")+0);",1}};
  
  for(const auto& [script,expected] : compoundChecks)
    {
      const auto tree=
        createParseTree(parser,script);
      const auto ast=
        executor.execParseTree(tree);
      observed=-1;
      std::visit(evaluator,*ast);
      assert(observed==expected);
    }
  
  // Registered host arithmetic supplies compound operations without native compound operators
  const std::pair<const char*,int> hostArithmeticChecks[]{
    {"a=ExampleNumber(20); b=a-3; record(b.value);",17},
    {"a=ExampleNumber(20); a-=3; record(a.value);",17},
    {"a=ExampleNumber(20); b=a/4; record(b.value);",5},
    {"a=ExampleNumber(20); a/=4; record(a.value);",5},
    {"a=ExampleNumber(20); a-=3; a/=2; record(a.value);",8}};
  
  for(const auto& [script,expected] : hostArithmeticChecks)
    {
      const auto tree=
        createParseTree(parser,script);
      const auto ast=
        executor.execParseTree(tree);
      observed=-1;
      std::visit(evaluator,*ast);
      assert(observed==expected);
    }
  
  // Reference arguments mutate caller cells while ordinary arguments remain values
  const std::pair<const char*,int> referenceChecks[]{
    {"a=3; fun f(&x){x=4;} f(&a); record(a);",4},
    {"a=3; fun f(x){x=4;} f(a); record(a);",3},
    {"a=3; fun f(&x){x=x+1;} fun relay(&y){f(&y);} relay(&a); record(a);",4},
    {"a=3; r=&a; fun f(&x){x=4;} f(r); record(a);",4},
    {"struct S{x;}; s=S(1); fun f(&x){x=7;} f(&s.x); record(s.x);",7},
    {"a=makeVec(1,2); fun f(&x){x=7;} f(&a[1]); record(a[1]);",7},
    {"fun init(&x){x=5;} init(&a); record(a);",5}};
  
  for(const auto& [script,expected] : referenceChecks)
    {
      const auto tree=
        createParseTree(parser,script);
      const auto ast=
        executor.execParseTree(tree);
      observed=-1;
      std::visit(evaluator,*ast);
      assert(observed==expected);
    }
  
  // Variadic packs distinguish copied values from shared references, including empty packs
  const std::pair<const char*,int> variadicChecks[]{
    {"fun f(...){return 1;} record(f());",1},
    {"fun f(...){return __VA__ARGS__[0]*10+__VA__ARGS__[1];} record(f(2,3));",23},
    {"fun f(x,...){return x*100+__VA__ARGS__[0]*10+__VA__ARGS__[1];} record(f(1,2,3));",123},
    {"fun f(x=7,...){return x;} record(f());",7},
    {"a=1; fun f(...){__VA__ARGS__[0]=9;} f(a); record(a);",1},
    {"fun f(&...){return 1;} record(f());",1},
    {"a=1; b=2; fun f(&...){__VA__ARGS__[0]+=2; __VA__ARGS__[1]*=3;} f(&a,&b); record(a*10+b);",36},
    {"a=1; fun f(x,&...){__VA__ARGS__[0]=x;} f(7,&a); record(a);",7},
    {"a=1; fun f(&x,...){x=__VA__ARGS__[0];} f(&a,5); record(a);",5},
    {"a=1; fun inner(&x){x=6;} fun outer(&...){inner(&__VA__ARGS__[0]);} outer(&a); record(a);",6},
    {"a=1; fun f(&...){return __VA__ARGS__;} refs=f(&a); refs[0]=9; record(a);",9},
    {"f=lambda(...){return __VA__ARGS__[0];}; record(f(5));",5}};
  
  for(const auto& [script,expected] : variadicChecks)
    {
      const auto tree=
        createParseTree(parser,script);
      const auto ast=
        executor.execParseTree(tree);
      observed=-1;
      std::visit(evaluator,*ast);
      assert(observed==expected);
    }
  
  // Valid integer arithmetic retains truncation, remainder signs, promotions, and floating-point behavior
  const std::pair<const char*,int> integerArithmeticChecks[]{
    {"record(7/3);",2},
    {"record(-7/3);",-2},
    {"record(-7%3);",-1},
    {"record(7%-3);",1},
    {"record(minimumInteger()/1);",std::numeric_limits<int>::min()},
    {"record(minimumInteger()%1);",0},
    {"x=minimumInteger(); x/=1; record(x);",std::numeric_limits<int>::min()},
    {"x=-7; x/=3; record(x);",-2},
    {"record((1/0.5==2.0)+0);",1},
    {"record((1.0/0>1.0)+0);",1},
    {"x=7; x/=2.0; record(x);",3},
    {"x=ExampleNumber(7); x.value/=2; record(x.value);",3}};
  
  for(const auto& [script,expected] : integerArithmeticChecks)
    {
      const auto tree=
        createParseTree(parser,script);
      const auto ast=
        executor.execParseTree(tree);
      observed=0;
      std::visit(evaluator,*ast);
      assert(observed==expected);
    }
  
  assert(DivisionOperator{}(static_cast<signed char>(-128),-1)==128);
  assert(RemainderOperator{}(static_cast<signed char>(-128),-1)==0);
  assert(DivisionOperator{}(std::numeric_limits<unsigned int>::max(),-1)==1u);
  
  // Named rejection modes verify division diagnostics in ordinary, compound, and member operations
  const std::pair<const char*,const char*> invalidArithmeticChecks[]{
    {"integer-division-zero","record(1/0);"},
    {"integer-remainder-zero","record(1%0);"},
    {"integer-division-overflow","record(minimumInteger()/-1);"},
    {"integer-remainder-overflow","record(minimumInteger()%-1);"},
    {"compound-division-zero","x=1; x/=0;"},
    {"compound-division-overflow","x=minimumInteger(); x/=-1;"},
    {"member-division-zero","x=ExampleNumber(1); x.value/=0;"},
    {"member-division-overflow","x=ExampleNumber(minimumInteger()); x.value/=-1;"}};
  
  for(const auto& [mode,script] : invalidArithmeticChecks)
    if(argc>1 and std::string(argv[1])==mode)
      {
        const auto tree=
          createParseTree(parser,script);
        const auto ast=
          executor.execParseTree(tree);
        std::visit(evaluator,*ast);
        errorEmitter("invalid integer arithmetic was unexpectedly accepted");
      }
  
  // Mixed positional/named binding, defaults, and variadic arguments remain valid
  const std::pair<const char*,int> argumentBindingChecks[]{
    {"fun f(x,y){return x*10+y;} record(f(1,.y=2));",12},
    {"fun f(x,y){return x*10+y;} record(f(.y=2,1));",12},
    {"fun f(x,y=2){return x*10+y;} record(f(.x=1));",12},
    {"fun f(x,...){return x*10+__VA__ARGS__[0];} record(f(1,2));",12},
    {"struct S{x;y;}; s=S(1,.y=2); record(s.x*10+s.y);",12},
    {"struct S{x;y;}; s=S(.y=2,1); record(s.x*10+s.y);",12},
    {"struct S{x;y=2;}; s=S(.x=1); record(s.x*10+s.y);",12}};
  
  for(const auto& [script,expected] : argumentBindingChecks)
    {
      const auto tree=
        createParseTree(parser,script);
      const auto ast=
        executor.execParseTree(tree);
      observed=-1;
      std::visit(evaluator,*ast);
      assert(observed==expected);
    }
  
  // The same parameter may not be filled twice, regardless of argument source order
  const std::pair<const char*,const char*> duplicateArgumentChecks[]{
    {"duplicate-function-mixed","fun f(x){return x;} f(1,.x=2);"},
    {"duplicate-function-reversed","fun f(x){return x;} f(.x=2,1);"},
    {"duplicate-function-named","fun f(x){return x;} f(.x=1,.x=2);"},
    {"duplicate-struct-mixed","struct S{x;}; S(1,.x=2);"},
    {"duplicate-struct-reversed","struct S{x;}; S(.x=2,1);"},
    {"duplicate-struct-named","struct S{x;}; S(.x=1,.x=2);"}};
  
  for(const auto& [mode,script] : duplicateArgumentChecks)
    if(argc>1 and std::string(argv[1])==mode)
      {
        const auto tree=
          createParseTree(parser,script);
        const auto ast=
          executor.execParseTree(tree);
        std::visit(evaluator,*ast);
        errorEmitter("duplicate argument was unexpectedly accepted");
      }
  
  if(argc>1 and std::string(argv[1])=="invalid-variadic-reference")
    {
      const auto tree=
        createParseTree(parser,"a=1; fun f(&...){return 1;} f(&a,2);");
      const auto ast=
        executor.execParseTree(tree);
      std::visit(evaluator,*ast);
      errorEmitter("non-reference variadic argument was unexpectedly accepted");
    }
  
  if(argc>1 and std::string(argv[1])=="uninitialized-compound")
    {
      const auto tree=
        createParseTree(parser,"x+=1;");
      const auto ast=
        executor.execParseTree(tree);
      std::visit(evaluator,*ast);
      errorEmitter("uninitialized compound assignment was unexpectedly accepted");
    }
  
  if(argc>1 and std::string(argv[1])=="invalid-compound-types")
    {
      const auto tree=
        createParseTree(parser,"x=1; x+=\"a\";");
      const auto ast=
        executor.execParseTree(tree);
      std::visit(evaluator,*ast);
      errorEmitter("invalid compound assignment was unexpectedly accepted");
    }
  
  if(argc>1 and std::string(argv[1])=="temporary-reference")
    {
      const auto tree=
        createParseTree(parser,"fun f(&x){x=4;} f(&(1+2));");
      const auto ast=
        executor.execParseTree(tree);
      std::visit(evaluator,*ast);
      errorEmitter("temporary reference was unexpectedly accepted");
    }
  
  if(argc>1 and std::string(argv[1])=="empty-function-result")
    {
      const auto tree=
        createParseTree(parser,"fun f(){x=1;} record(f());");
      const auto ast=
        executor.execParseTree(tree);
      std::visit(evaluator,*ast);
      errorEmitter("empty function result was unexpectedly accepted");
    }
}

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

/// Check launch ordering, literal arguments, diagnostics, and ownership before file reading
void validateLaunchArgs(const int narg,
                        const char* const* arg)
{
  // Invalid command lines are checked in separate processes for message and status
  if(narg>1)
    {
      const std::string_view mode=arg[1];
      
      if(mode=="launch-missing-file")
        {
          const char* invalid[]{"nice","setup();","-f"};
          (void)parseLaunchArgs(3,invalid);
        }
      else if(mode=="launch-unknown-option")
        {
          const char* invalid[]{"nice","-e","setup();"};
          (void)parseLaunchArgs(3,invalid);
        }
      else if(mode=="launch-attached-path")
        {
          const char* invalid[]{"nice","-fcommon.nice"};
          (void)parseLaunchArgs(2,invalid);
        }
    }
  
  const char* mixed[]{"nice","setup();","-f","common.nice","finish();","-f","last.nice"};
  const auto parts=parseLaunchArgs(7,mixed);
  assert(parts.size()==4);
  assert(parts[0].kind==ProgramPart::Kind::Immediate and parts[0].value=="setup();" and parts[0].argIndex==1);
  assert(parts[1].kind==ProgramPart::Kind::File and parts[1].value=="common.nice" and parts[1].argIndex==3);
  assert(parts[2].kind==ProgramPart::Kind::Immediate and parts[2].value=="finish();" and parts[2].argIndex==4);
  assert(parts[3].kind==ProgramPart::Kind::File and parts[3].value=="last.nice" and parts[3].argIndex==6);
  
  // After --, every argument is program text, including option-looking text
  const char* literal[]{"nice","--","-f","--","-1;",""};
  const auto literalParts=parseLaunchArgs(6,literal);
  assert(literalParts.size()==4);
  
  for(std::size_t i=0;i<literalParts.size();i++)
    {
      assert(literalParts[i].kind==ProgramPart::Kind::Immediate);
      assert(literalParts[i].value==literal[i+2] and literalParts[i].argIndex==i+2);
    }
  
  // A file operand is taken literally, while empty launches remain valid
  const char* path[]{"nice","-f","--","record();"};
  const auto pathParts=parseLaunchArgs(4,path);
  assert(pathParts.size()==2 and pathParts[0].kind==ProgramPart::Kind::File);
  assert(pathParts[0].value=="--" and pathParts[0].argIndex==2);
  const char* empty[]{"nice"};
  assert(parseLaunchArgs(1,empty).empty());
  
  // Parsed descriptors own their values independently of arg storage
  char text[]{"value=1;"};
  const char* owned[]{"nice",text};
  const auto ownedParts=parseLaunchArgs(2,owned);
  text[0]='X';
  assert(ownedParts[0].value=="value=1;");
  
}

/// Check master-only file reading, owned text, shared scope, and compilation equivalence
void validateProgramPieces(const int narg,
                           const char* const* arg)
{
  if(narg>1)
    {
      const std::string_view mode=arg[1];
      
      if(mode=="program-missing-file")
        (void)readAndBroadcastProgramPieces({{ProgramPart::Kind::File,"/dev/null/missing.nice",3}});
      else if(mode=="program-unreadable-file")
        (void)readAndBroadcastProgramPieces({{ProgramPart::Kind::File,"/",3}});
    }
  
  // Library broadcasts preserve full-width counts and binary strings from any root
  const int root=nissa::nranks-1;
  const uint64_t expected=std::numeric_limits<uint64_t>::max()-3;
  const uint64_t local=nissa::rank==root?expected:0;
  assert(nissa::broadcast(local,root)==expected);
  assert(nissa::broadcast(nissa::rank==root?17:0,root)==17);
  std::string binary=nissa::rank==root?std::string("a\0bc",4):"stale receiver text";
  nissa::broadcast(binary,root);
  assert(binary==std::string("a\0bc",4));
  std::string empty=nissa::rank==root?"":"stale receiver text";
  nissa::broadcast(empty,root);
  assert(empty.empty());
  
  // The default source follows master_rank, including when it is not rank zero
  const int savedMaster=nissa::master_rank;
  nissa::master_rank=root;
  assert(nissa::broadcast(local)==expected);
  assert(nissa::broadcast(nissa::rank==root?2.5:0.0)==2.5);
  assert(nissa::broadcast(nissa::rank==root?1.25f:0.0f)==1.25f);
  std::string defaultRoot=nissa::rank==root?"master":"other";
  nissa::broadcast(defaultRoot);
  assert(defaultRoot=="master");
  nissa::master_rank=savedMaster;
  
  const std::string prefix="corrBlock(\"pion\"); value=1;";
  const std::string fileText=std::string(70000,' ')+
    "value+=2; { hidden=99; value+=3; } fun increment(){return value+1;} line=Gamma(5)*Source(0);";
  const std::string suffix=
    "record(increment()); corrGammas(\"pion\",5,5); corrTrace(\"pion\",line,line);";
  LaunchArgs parts;
  char path[]{"/tmp/nice-program-XXXXXX"};
  
  if(nissa::rank==nissa::master_rank)
    {
      const int descriptor=mkstemp(path);
      assert(descriptor>=0);
      FILE* file=fdopen(descriptor,"wb");
      assert(file!=nullptr);
      const auto written=std::fwrite(fileText.data(),1,fileText.size(),file);
      assert(written==fileText.size());
      const int status=std::fclose(file);
      assert(status==0);
      parts={{ProgramPart::Kind::Immediate,prefix,1},
             {ProgramPart::Kind::File,path,3},
             {ProgramPart::Kind::File,"/dev/null",5},
             {ProgramPart::Kind::Immediate,"",6},
             {ProgramPart::Kind::Immediate,suffix,7}};
    }
  else
    // Non-master descriptors are deliberately unusable: only master supplies sources
    parts={{ProgramPart::Kind::File,"/dev/null/non-master-must-not-open",99}};
  
  const auto pieces=readAndBroadcastProgramPieces(parts);
  
  if(nissa::rank==nissa::master_rank)
    {
      const int status=std::remove(path);
      assert(status==0);
    }
  
  assert(pieces.size()==5);
  assert(pieces[0].text==prefix and pieces[0].origin=="arg[1]");
  assert(pieces[1].text==fileText and pieces[1].origin.starts_with("arg[3] file '/tmp/nice-program-"));
  assert(pieces[2].text.empty() and pieces[2].origin=="arg[5] file '/dev/null'");
  assert(pieces[3].text.empty() and pieces[3].origin=="arg[6]");
  assert(pieces[4].text==suffix and pieces[4].origin=="arg[7]");
  assert(readAndBroadcastProgramPieces({}).empty());
  
  const auto evaluate=[](const std::vector<ProgramPiece>& sources)
  {
    // Parsing returns before evaluation, so the merged AST must retain function bodies
    const ASTNode ast=parseProgramPieces(sources);
    const auto& body=std::get<ASTNodesNode>(ast);
    assert(std::holds_alternative<ASTNodesNode>(*body.list[3]));
    Nice nice;
    Evaluator evaluator;
    setupNice(evaluator,nice.correlations);
    int observed{};
    registerHostFunction(evaluator,"record",{"value"},[&observed](const int value)
    {
      observed=value;
    });
    std::visit(evaluator,ast);
    assert(observed==7);
    
    return writeLegacyInput(nice.compile());
  };
  
  // The same requests supplied as files or immediate text produce identical instructions
  const std::string mixedOutput=evaluate(pieces);
  const auto immediate=readAndBroadcastProgramPieces({
    {ProgramPart::Kind::Immediate,prefix+fileText+suffix,1}});
  assert(evaluate(immediate)==mixedOutput);
  assert(parseProgramPieces({}).list.empty());
}

/// Run positive checks normally; a named argument selects an expected fatal diagnostic
int main(int narg,
         char** arg)
{
  nissa::init_MPI_thread(narg,arg);
  nissa::get_MPI_nranks();
  nissa::get_MPI_rank();
  
  validateLaunchArgs(narg,arg);
  validateProgramPieces(narg,arg);
  validateCorrParts();
  validateGammaAdjoints();
  validateTaskGraph(narg,arg);
  validateGraphCompiler(narg,arg);
  validateScheduling(narg,arg);
  validateLegacyWriter(narg,arg);
  validateScriptBindings(narg,arg);
  validateProgramGrammar(narg,arg);
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
  
#ifdef USE_MPI
  MPI_Finalize();
#endif
}
