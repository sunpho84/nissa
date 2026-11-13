#include "corrParts.hpp"
#include "configurationSources.hpp"
#include "correlations.hpp"
#include "correlationTaskGraph.hpp"
#include "graphCompiler.hpp"
#include "gaugeRecipes.hpp"
#include "lineValue.hpp"
#include "taskExecution.hpp"
#include "programExecution.hpp"
#include "schedulingStrategies.hpp"
#include "legacyInputWriter.hpp"
#include "niceSetup.hpp"
#include "nice.hpp"
#include "runOptions.hpp"
#include "nASTy.hpp"
#include "new_types/dirac.hpp"
#include "base/multiGridParams.hpp"
#include "base/DDalphaAMG_bridge.hpp"
#include "dirac_operators/tmclovQ/dirac_operator_tmclovQ.hpp"
#include "operations/su3_paths/clover_term.hpp"
#include "routines/ios.hpp"
#include "routines/mpi_routines.hpp"
#include "base/init.hpp"
#include "base/close.hpp"
#include "io/writer.hpp"
#include "linalgs/reduce.hpp"
#include "operations/gaugeconf.hpp"
#include "operations/su3_paths/plaquette.hpp"

#include <algorithm>
#include <cassert>
#include <complex>
#include <filesystem>
#include <format>
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
    "CloverRun 0\n\n"
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
  
  Correlations cloverRequests;
  cloverRequests.addBlock("clover").addTrace(
    Prop{.r=1,.csw=1.25}*source,Prop{.mass=0.2,.r=-1,.csw=1.25}*source);
  auto cloverCompiled=scheduleCorrelations(buildTaskGraph(cloverRequests));
  assert(writeLegacyInput(cloverCompiled).starts_with("CloverRun 1\ncSW 1.25\n\n"));
  
  if(argc>1 and std::string_view(argv[1])=="legacy-mixed-csw")
    {
      std::get<Prop>(std::get<ExtendSelectingTPars>(cloverCompiled.instructions[2].operation).pars).csw=0.5;
      writeLegacyInput(cloverCompiled);
      errorEmitter("mixed clover coefficients were unexpectedly serialized");
    }
  
  if(argc>1 and std::string_view(argv[1])=="legacy-nonfinite-csw")
    {
      std::get<Prop>(std::get<ExtendSelectingTPars>(cloverCompiled.instructions[1].operation).pars).csw=
        std::numeric_limits<double>::infinity();
      writeLegacyInput(cloverCompiled);
      errorEmitter("nonfinite clover coefficient was unexpectedly serialized");
    }
  
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
  
  const auto cloverTree=createParseTree(parser,R"(
    setDefaultPropCsw(1.75);
    firstClover=Prop(.r=1);
    recordProp(firstClover);
    recordProp(Prop(.r=1,.csw=0.5));
    setDefaultPropCsw(0.0);
    recordProp(firstClover);
  )");
  const auto cloverAst=defaultsExecutor.execParseTree(cloverTree);
  std::visit(evaluator,*cloverAst);
  assert(observedProps.size()==12);
  assert(observedProps[9].csw==1.75 and observedProps[10].csw==0.5);
  assert(observedProps[11].csw==1.75 and Prop{}.csw==0.0);
  assert(observedProps[9].dag().csw==1.75);
  
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
      else if(mode=="program-missing-run-settings")
        (void)prepareProgram({{ProgramPart::Kind::Immediate,"",1}});
      else if(mode=="program-unreadable-file")
        (void)readAndBroadcastProgramPieces({{ProgramPart::Kind::File,"/",3}});
    }
  
  // Library broadcasts preserve full-width counts and binary strings from any root
  const int root=nissa::nranks-1;
  const uint64_t expected=std::numeric_limits<uint64_t>::max()-3;
  const uint64_t local=nissa::rank==root?expected:0;
  assert(nissa::broadcast(local,root)==expected);
  assert(nissa::broadcast(nissa::rank==root?17:0,root)==17);
  assert(nissa::broadcast(nissa::rank==root,root));
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

  // Read a file larger than one chunk, with an unusable worker path in collective modes
  const std::string readPath=nissa::is_master_rank()?path:"/dev/null/non-master-must-not-open";
  const std::string masterBytes=nissa::readFileBytes(readPath,nissa::FileReadMode::MasterOnly);
  if(nissa::is_master_rank())
    assert(masterBytes==fileText);
  else
    assert(masterBytes.empty());
  
  const std::string sharedBytes=nissa::readFileBytes(readPath,nissa::FileReadMode::Collective);
  assert(sharedBytes==fileText);
  
  const std::string independentPath=nissa::is_master_rank()?path:"/dev/null";
  const std::string independentBytes=
    nissa::readFileBytes(independentPath,nissa::FileReadMode::Independent);
  if(nissa::is_master_rank())
    assert(independentBytes==fileText);
  else
    assert(independentBytes.empty());
  
  const auto pieces=readAndBroadcastProgramPieces(parts);
  
  // Prepare settings and correlations together while the file still exists
  LaunchArgs runParts=parts;
  runParts.insert(runParts.begin(),{ProgramPart::Kind::Immediate,
    "fun record(value){return value;} setGeometry(24,48); setRun(7,8,2); setDilution(1,0);",8});
  const PreparedProgram prepared=prepareProgram(runParts);
  
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
  
  // Returned instructions remain usable after preparation's AST and evaluator are destroyed
  prepared.options.validate();
  prepared.correlations.validate();
  assert(prepared.options.L==24 and prepared.options.T==48);
  assert(prepared.options.seed==7 and prepared.options.nHits==8 and prepared.options.nCopies==2);
  assert(prepared.options.diluteSpin==true and prepared.options.diluteColor==false);
  assert(writeLegacyInput(prepared.correlations)==mixedOutput);
  const auto immediate=readAndBroadcastProgramPieces({
    {ProgramPart::Kind::Immediate,prefix+fileText+suffix,1}});
  assert(evaluate(immediate)==mixedOutput);
  assert(parseProgramPieces({}).list.empty());
}

/// Check list rows, path resolution, ordered source registration, and owned preparation results
void validateConfigurationSources(const int narg,
                                  const char* const* arg)
{
  if(narg>1)
    {
      const std::string_view mode=arg[1];
      
      if(mode=="list-malformed-row")
        (void)parseConfigurationList("conf output extra\n","/tmp/configurations.list");
      else if(mode=="list-extra-column")
        (void)parseConfigurationList("conf\toutput\textra\n","/tmp/configurations.list");
      else if(mode=="list-empty-tab-field")
        (void)parseConfigurationList("conf\t   \n","/tmp/configurations.list");
      else if(mode=="list-missing-file")
        (void)materializeAndBroadcastConfigurations({{"/dev/null/missing.list"}});
      else if(mode=="list-duplicate-output")
        {
          char firstPath[]{"/tmp/nice-duplicate-first-XXXXXX"};
          char secondPath[]{"/tmp/nice-duplicate-second-XXXXXX"};
          const int firstDescriptor=mkstemp(firstPath);
          const int secondDescriptor=mkstemp(secondPath);
          assert(firstDescriptor>=0 and secondDescriptor>=0);
          FILE* firstFile=fdopen(firstDescriptor,"wb");
          FILE* secondFile=fdopen(secondDescriptor,"wb");
          assert(firstFile!=nullptr and secondFile!=nullptr);
          assert(std::fputs("a\tout/../same\n",firstFile)>=0);
          assert(std::fputs("b\tsame\n",secondFile)>=0);
          assert(std::fclose(firstFile)==0 and std::fclose(secondFile)==0);
          (void)materializeAndBroadcastConfigurations({{firstPath,secondPath}});
        }
    }
  
  char path[]{"/tmp/nice-list-XXXXXX"};
  char secondPath[]{"/tmp/nice-list-next-XXXXXX"};
  
  if(nissa::rank==nissa::master_rank)
    {
      const int descriptor=mkstemp(path);
      assert(descriptor>=0);
      FILE* file=fdopen(descriptor,"wb");
      assert(file!=nullptr);
      assert(std::fputs("  # comment\r\n\r\ninput one\toutput one\r\nsub/input   results/two\n",file)>=0);
      assert(std::fclose(file)==0);
      const int secondDescriptor=mkstemp(secondPath);
      assert(secondDescriptor>=0);
      FILE* secondFile=fdopen(secondDescriptor,"wb");
      assert(secondFile!=nullptr);
      assert(std::fputs("third\toutput three\nfourth output_four\n",secondFile)>=0);
      assert(std::fclose(secondFile)==0);
    }
  
  // Only the master owns the list path; the program itself distributes that path
  const std::string script=nissa::rank==nissa::master_rank?
    "setGeometry(16,32); configurationsFromFile(\""+std::string(path)+"\"); "
    "configurationsFromFile(\""+std::string(secondPath)+"\");":"invalid local script";
  const PreparedProgram prepared=prepareProgram({{ProgramPart::Kind::Immediate,script,1}});
  const ConfigurationSources masterOnly=nissa::rank==nissa::master_rank?
    ConfigurationSources{{path}}:ConfigurationSources{{"/dev/null/non-master-must-not-open"}};
  const auto masterEntries=materializeAndBroadcastConfigurations(masterOnly);
  assert(masterEntries.size()==2);
  std::string firstSource=nissa::rank==nissa::master_rank?path:"";
  std::string secondSource=nissa::rank==nissa::master_rank?secondPath:"";
  nissa::broadcast(firstSource);
  nissa::broadcast(secondSource);
  
  if(nissa::rank==nissa::master_rank)
    {
      assert(std::remove(path)==0);
      assert(std::remove(secondPath)==0);
    }
  
  assert(prepared.configurationEntries.size()==4);
  const std::string workingDirectory=std::filesystem::current_path().string();
  assert(prepared.configurationEntries[0].path==workingDirectory+"/input one");
  assert(prepared.configurationEntries[0].outputDirectory==workingDirectory+"/output one");
  assert(prepared.configurationEntries[0].sourceList==firstSource);
  assert(prepared.configurationEntries[0].sourceRow==3);
  assert(prepared.configurationEntries[1].path==workingDirectory+"/sub/input");
  assert(prepared.configurationEntries[1].outputDirectory==workingDirectory+"/results/two");
  assert(prepared.configurationEntries[1].sourceList==firstSource);
  assert(prepared.configurationEntries[1].sourceRow==4);
  assert(prepared.configurationEntries[2].path==workingDirectory+"/third");
  assert(prepared.configurationEntries[2].outputDirectory==workingDirectory+"/output three");
  assert(prepared.configurationEntries[2].sourceList==secondSource);
  assert(prepared.configurationEntries[2].sourceRow==1);
  assert(prepared.configurationEntries[3].path==workingDirectory+"/fourth");
  assert(prepared.configurationEntries[3].outputDirectory==workingDirectory+"/output_four");
  assert(prepared.configurationEntries[3].sourceList==secondSource);
  assert(prepared.configurationEntries[3].sourceRow==2);
}

/// Check deferred run settings, complete-group overrides, and script argument validation
void validateRunSettings(const int narg,
                         const char* const* arg)
{
  RunOptions options;
  Evaluator evaluator;
  registerRunOptions(evaluator,options);
  const auto execute=[&](const std::vector<ProgramPiece>& pieces)
  {
    const ASTNode ast=parseProgramPieces(pieces);
    std::visit(evaluator,ast);
  };
  
  // Direct C++ callers use the same data and validation without constructing a script
  const RunOptions defaults{.L=16,.T=32};
  defaults.validate();
  assert(defaults.seed==Sitmo::Rng::DEFAULT_SEED());
  assert(defaults.diluteSpin==false and defaults.diluteColor==false);
  RunOptions largeGeometry=defaults;
  largeGeometry.L=200;
  largeGeometry.T=200;
  largeGeometry.validate();
  assert(defaults.nHits==1 and defaults.nCopies==1);
  assert(defaults.temporalBoundary==RunOptions::TemporalBoundary::Antiperiodic);
  
  if(narg>1)
    {
      const std::string_view mode=arg[1];
      const auto reject=[&](const std::string_view selected,
                           const auto& change)
      {
        if(mode==selected)
          {
            RunOptions invalid=defaults;
            change(invalid);
            invalid.validate();
            CRASH("Invalid run settings were unexpectedly accepted");
          }
      };
      reject("run-invalid-geometry",[](RunOptions& value)
      {
        value.L=0;
      });
      reject("run-overflow-volume",[](RunOptions& value)
      {
        value.L=std::numeric_limits<int>::max();
        value.T=std::numeric_limits<int>::max();
      });
      reject("run-signed-volume",[](RunOptions& value)
      {
        value.L=2097152;
        value.T=1;
      });
      reject("run-invalid-seed",[](RunOptions& value)
      {
        value.seed=-1;
      });
      reject("run-oversized-seed",[](RunOptions& value)
      {
        value.seed=uint64_t{std::numeric_limits<uint32_t>::max()}+1;
      });
      reject("run-invalid-hits",[](RunOptions& value)
      {
        value.nHits=0;
      });
      reject("run-invalid-copies",[](RunOptions& value)
      {
        value.nCopies=0;
      });
      reject("run-negative-origin",[](RunOptions& value)
      {
        value.drawGlobalOriginTask.enabled=false;
        value.drawGlobalOriginTask.fixedOrigin[1]=-1;
      });
      reject("run-origin-time-outside",[](RunOptions& value)
      {
        value.drawGlobalOriginTask.enabled=false;
        value.drawGlobalOriginTask.fixedOrigin[0]=value.T;
      });
      reject("run-origin-space-outside",[](RunOptions& value)
      {
        value.drawGlobalOriginTask.enabled=false;
        value.drawGlobalOriginTask.fixedOrigin[3]=value.L;
      });
      reject("run-invalid-boundary-value",[](RunOptions& value)
      {
        value.temporalBoundary=static_cast<RunOptions::TemporalBoundary>(17);
      });
      
      if(mode=="run-invalid-dilution")
        execute({{"setDilution(2,0);","invalid dilution"}});
      else if(mode=="run-invalid-boundary")
        execute({{"setTemporalBoundary(\"other\");","invalid boundary"}});
      else if(mode=="rng-negative-script-limit")
        execute({{"reserveRngDoublesPerSite(-1);","invalid RNG limit"}});
      else if(mode=="run-invalid-origin")
        execute({{"setRandomOriginEnabled(2);","invalid random origin setting"}});
      else if(mode=="run-negative-script-seed" or mode=="run-oversized-script-seed")
        {
          execute({{mode=="run-negative-script-seed"?
                     "setGeometry(4,8); setRun(-1,1,1);":
                     "setGeometry(4,8); setRun(4294967296,1,1);","invalid script seed"}});
          options.validate();
        }
      else if(mode=="script-integer-literal-overflow")
        execute({{"value=9223372036854775808;","oversized integer literal"}});
      else if(mode=="script-integer-literal-huge")
        execute({{"value=999999999999999999999999999999999999;","oversized integer literal"}});
      else if(mode=="run-wide-hit-count")
        execute({{"setRun(7,2147483648,1);","oversized host int argument"}});
    }
  
  // Geometry alone inherits Sitmo's seed and no source dilution
  execute({{"setGeometry(16,32);","default settings"}});
  options.validate();
  assert(options.seed==Sitmo::Rng::DEFAULT_SEED());
  assert(options.diluteSpin==false and options.diluteColor==false);
  
  // Both integer representations reach the full unsigned Sitmo seed range without truncation
  for(const int64_t seed : {int64_t{0},int64_t{2147483647},int64_t{2147483648},
                            int64_t{3472291050},int64_t{4294967295}})
    {
      execute({{std::format("setRun(.seed={},.nHits=1,.nCopies=1);",seed),"full-range seed"}});
      options.validate();
      assert(options.seed==seed);
    }
  
  // Large literals retain their exact value through ordinary script storage and formatting
  std::vector<int64_t> wideIntegers;
  std::string wideText;
  registerHostFunction(evaluator,"recordWideInteger",{"value"},[&wideIntegers](const int64_t value)
  {
    wideIntegers.push_back(value);
  });
  registerHostFunction(evaluator,"recordWideText",{"value"},[&wideText](const std::string& value)
  {
    wideText=value;
  });
  execute({{"largeInteger=9223372036854775807; recordWideInteger(largeInteger); "
            "recordWideText(to_string(largeInteger)); recordWideInteger(4294967295+1); "
            "recordWideInteger(seq(2)[2147483648-2147483647]);",
            "wide script integers"}});
  assert((wideIntegers==std::vector<int64_t>{std::numeric_limits<int64_t>::max(),4294967296,1}));
  assert(wideText=="9223372036854775807");
  
  // Earlier incomplete settings can be replaced before validation, with no setup side effects
  execute({{"setGeometry(-1,-2); setRun(-1,0,0); setDilution(0,0);","file defaults"},
           {"setGeometry(.T=48,.L=24); setRun(.nCopies=2,.seed=1234,.nHits=8); "
            "setDilution(.diluteColor=0,.diluteSpin=1); setTemporalBoundary(\"periodic\"); "
            "setRandomOriginEnabled(0);","argument override"}});
  options.validate();
  assert(options.L==24 and options.T==48 and options.seed==1234);
  assert(options.nHits==8 and options.nCopies==2);
  assert(not options.drawGlobalOriginTask.enabled);
  assert(options.reservedRngDoublesPerSite==1'000'000);
  assert(options.diluteSpin==true and options.diluteColor==false);
  assert(options.temporalBoundary==RunOptions::TemporalBoundary::Periodic);
  
  // Each repeated setter replaces the full group and keeps unrelated settings
  execute({{"setRun(7,1,3); setGeometry(12,24); setDilution(0,1); "
            "setRandomOriginEnabled(1); "
            "setTemporalBoundary(\"antiperiodic\");","later override"}});
  options.validate();
  assert(options.L==12 and options.T==24 and options.seed==7);
  assert(options.nHits==1 and options.nCopies==3);
  assert(options.drawGlobalOriginTask.enabled);
  assert(options.reservedRngDoublesPerSite==1'000'000);
  assert(options.diluteSpin==false and options.diluteColor==true);
  assert(options.temporalBoundary==RunOptions::TemporalBoundary::Antiperiodic);
  
  execute({{"setRun(7,8,2); reserveRngDoublesPerSite(.reservedDoublesPerSite=12);", "RNG reservation"}});
  options.validate();
  assert(options.reservedRngDoublesPerSite==12);
  
  // Fixed coordinates are checked against the final geometry and consume no draws
  execute({{"setOrigin(.z=1,.t=7,.y=2,.x=3); setGeometry(4,8);","fixed origin"}});
  options.validate();
  const nissa::Coords fixedOrigin{7,3,2,1};
  assert(not options.drawGlobalOriginTask.enabled);
  assert(options.drawGlobalOriginTask.fixedOrigin==fixedOrigin);
  assert(options.drawGlobalOriginTask.doubleDrawCost()==0);
  
  // Random selection can be reenabled without discarding the explicitly selected coordinates
  execute({{"setRandomOriginEnabled(1);","random origin"}});
  assert(options.drawGlobalOriginTask.enabled);
  assert(options.drawGlobalOriginTask.doubleDrawCost()==NDIM);
  execute({{"setRandomOriginEnabled(0);","restore fixed origin"}});
  options.validate();
  assert(options.drawGlobalOriginTask.fixedOrigin==fixedOrigin);
}

/// Check deferred recipe construction, empty storage, and invalid recipe diagnostics
void validateGaugeRecipeBindings(const int narg,
                                 const char* const* arg)
{
  if(narg>1 and std::string_view(arg[1])=="gauge-unavailable-scheduled-recipe")
    {
      (void)prepareProgram({{ProgramPart::Kind::Immediate,
        "setGeometry(4,4); raw=loadedConf(); corrBlock(\"invalid\"); "
        "corrGammas(\"invalid\",5,5); source=Source(); "
        "line=Prop(.gaugeRecipe=1,.r=1)*source; corrTrace(\"invalid\",line,line);",1}});
      errorEmitter("unavailable scheduled gauge recipe was unexpectedly accepted");
    }
  
  PreparedProgram prepared=prepareProgram({{ProgramPart::Kind::Immediate,
    "setGeometry(4,4); rawRecipe=loadedConf(); copy=rawRecipe; "
    "zero=apeSmearConf(.inputRecipe=copy,.alpha=0.4,.n=0); "
    "smooth=apeSmearConf(rawRecipe,0.4,2); "
    "chained=apeSmearConf(smooth,0.2,1);",1}});
  auto& recipes=prepared.gaugeRecipes;
  assert(not nissa::lxGeomInited);
  assert(recipes.size()==4);
  
  // Script evaluation records producers and reserves storage without allocating gauge fields
  for(const auto& recipe : recipes)
    assert(not recipe.field);
  
  if(narg>1)
    {
      const std::string_view mode=arg[1];
      
      if(mode=="gauge-invalid-handle")
        apeSmearConf(recipes,recipes.size(),0.4,1);
      else if(mode=="gauge-negative-steps")
        apeSmearConf(recipes,0,0.4,-1);
      else if(mode=="gauge-nonfinite-alpha")
        apeSmearConf(recipes,0,std::numeric_limits<double>::infinity(),1);
      else if(mode=="gauge-empty-producer")
        {
          recipes[1].produce={};
          validateGaugeRecipes(recipes);
        }
      else if(mode=="gauge-unavailable-field")
        (void)gaugeField(recipes,0);
      else if(mode=="gauge-before-geometry")
        prepareGaugeFields(recipes,{},Sitmo::Rng(1),
                           RngLayout{prepared.options},0);
    }
}

/// Check runtime dispatch when an output-first schedule reorders independent tasks
void validateTaskExecution(const RunOptions& options)
{
  const Source first(0);
  const Source second(0);
  const std::set<std::pair<int,int>> gammas{{5,5}};
  const CorrelationTaskGraph graph(
    {Task{.id=0,.operation=first,.inputs{},.labels{}},
     Task{.id=1,.operation=second,.inputs{},.labels{}},
     Task{.id=2,.operation=ExtendSelectingTPars{Pars(Gamma{.iGamma=1})},
          .inputs{{0,2.0}},.labels{}},
     Task{.id=3,.operation=TraceContraction{gammas},
          .inputs{{2,1.0},{1,1.0}},.labels{}}},
    {OutputBlock{.name="dispatch",.gammaList=gammas,.traces{3}}});
  const CompiledCorrelations compiled=scheduleCorrelations(graph,scheduleByOutput);
  std::string callOrder;
  CompRegistry registry;
  const CompKind gammaKind=registry.intern("GammaPair");
  const CompKind timeKind=registry.intern("TimeSeparation");
  const TaskHandlers handlers{
    .source=[&](const Source&)
    {
      callOrder+='S';
      LineValue result(options,nissa::MemorySpace::CPU,nissa::SpaceTimeLayout::CPU);
      result.field.tensor.data[0]=static_cast<double>(callOrder.size());
      
      return result;
    },
    .extend=[&](const ExtendSelectingTPars&,
                const std::span<const WeightedLineInput> inputs)
    {
      callOrder+='E';
      assert(inputs.size()==1 and inputs[0].weight==2.0);
      LineValue result(options,nissa::MemorySpace::CPU,nissa::SpaceTimeLayout::CPU);
      result.field.tensor.data[0]=inputs[0].weight*inputs[0].value->field.tensor.data[0];
      
      return result;
    },
    .trace=[&](const TraceContraction& contraction,
               const LineValue& backward,
               const LineValue& forward)
    {
      callOrder+='T';
      assert(contraction.gammaList==gammas);
      Tensor result{{ScalarKind::Complex,
                     {{gammaKind,CompSide::None,contraction.gammaList.size()},
                      {StandardComp::Copy,CompSide::None,backward.nCopies},
                      {timeKind,CompSide::None,static_cast<size_t>(options.T)}}},
                    nissa::MemorySpace::CPU};
      result.data[0]=backward.field.tensor.data[0]+forward.field.tensor.data[0];
      
      return result;
    }};
  const TaskValues values=executeScheduledTasks(compiled,handlers);
  assert(callOrder=="SEST");
  assert(values.size()==4);
  assert(not values[0] and not values[1] and not values[2]);
  assert(std::get<Tensor>(*values[3]).data[0]==5.0);
}

/// Check volume-scaled budgets, independent scopes, direct skipping, and actual integer draws
void validateRngLayout(const int narg,
                       const char* const* arg)
{
  const auto equalPositions=
    [](const RngPosition& left,
       const RngPosition& right)
  {
    return std::equal(left.data,left.data+4,right.data);
  };
  
  // Each configuration and hit receives the same budget, scaled to a 200^4 lattice
  const RunOptions options{.L=200,.T=200};
  RngLayout layout(options);
  layout.origin=
    layout.addOperation(RngScope::Hit,GlobalOriginOperation{}.doubleDrawCost());
  const auto gauge=
    layout.addOperation(RngScope::GaugeConf,1'600'000'000'000);
  const auto deterministic=
    layout.addOperation(RngScope::GaugeConf,0);
  const auto secondGauge=
    layout.addOperation(RngScope::GaugeConf,17);
  const auto source=
    layout.addOperation(RngScope::Hit,32);
  const Sitmo::Rng rng(3472291050U);
  const auto origin=
    layout.operationRange(rng,layout.origin,3,5);
  const auto gaugeRange=
    layout.operationRange(rng,gauge,3);
  const auto secondRange=
    layout.operationRange(rng,secondGauge,3);
  const auto sourceRange=
    layout.operationRange(rng,source,3,5);
  
  // Gauge configuration preparation occupies slot zero; hit five occupies slot six
  assert(equalPositions(origin.position,RngPosition{{{9'600'000'000'000'000,4'800'000'000'000'000,0,0}}}));
  assert(equalPositions(gaugeRange.position,RngPosition{{{0,4'800'000'000'000'000,0,0}}}));
  assert(equalPositions(secondRange.position,RngPosition{{{1'600'000'000'000,4'800'000'000'000'000,0,0}}}));
  assert(equalPositions(sourceRange.position,origin.position+uint64_t{NDIM}));
  assert(layout.operationRange(rng,deterministic,3).doubleDrawCost==0);
  
  // An omitted or explicitly absent hit selects gauge preparation; hit zero remains a real hit
  assert(equalPositions(layout.operationRange(rng,gauge,3,std::nullopt).position,gaugeRange.position));
  
  // Compare doubles with independently addressed raw outputs, including both adjacent operations
  for(uint64_t i=0;i<17;i++)
    {
      const RngPosition raw{{{3'200'000'000'000+2*i,9'600'000'000'000'000,0,0}}};
      const uint32_t low=rng.generateNth(raw);
      const uint32_t high=rng.generateNth(raw+uint64_t{1});
      constexpr double F=4294967296.0;
      const double expected=std::min((low/F+high)/F,0.999999999999999889);
      assert(secondRange.uniformAt(i)==expected);
    }
  
  // Every hit receives its own origin, independently of execution order and requested hit count
  const auto nextHit=
    layout.operationRange(rng,layout.origin,3,6);
  assert(equalPositions(nextHit.position,origin.position+layout.reservedDoubles()));
  assert(equalPositions(layout.operationRange(rng,layout.origin,4,0).position,
                        RngPosition{{{1'600'000'000'000'000,6'400'000'000'000'000,0,0}}}));
  const auto skipped=
    layout.operationRange(rng,layout.origin,10000,7);
  const double retained=origin.uniformAt(0);
  assert(skipped.uniformAt(0)==layout.operationRange(rng,layout.origin,10000,7).uniformAt(0));
  assert(origin.uniformAt(0)==retained);
  assert(equalPositions(rng.state,RngPosition{}));
  
  // Deterministic operations consume no space, and configuration costs do not move hit operations
  RngLayout reordered(options);
  reordered.addOperation(RngScope::Hit,0);
  reordered.origin=
    reordered.addOperation(RngScope::Hit,NDIM);
  reordered.addOperation(RngScope::GaugeConf,1);
  assert(equalPositions(reordered.operationRange(rng,reordered.origin,3,5).position,origin.position));
  
  // A larger lattice scales all reserved ranges automatically
  RunOptions largerOptions=options;
  largerOptions.T*=2;
  RngLayout larger(largerOptions);
  larger.origin=
    larger.addOperation(RngScope::Hit,NDIM);
  assert(equalPositions(larger.operationRange(rng,larger.origin,3,5).position,origin.position*uint64_t{2}));
  
  // Both scopes may use their complete budget without sharing any double
  const RunOptions exactOptions{.L=1,.T=2,.reservedRngDoublesPerSite=2};
  RngLayout exact(exactOptions);
  const auto exactGauge=
    exact.addOperation(RngScope::GaugeConf,4);
  const auto exactHit=
    exact.addOperation(RngScope::Hit,4);
  assert(equalPositions(exact.operationRange(rng,exactHit,0,0).position,
                        exact.operationRange(rng,exactGauge,0).position+uint64_t{4}));
  assert(equalPositions(exact.operationRange(rng,exactHit,0,1).position,
                        exact.operationRange(rng,exactHit,0,0).position+uint64_t{4}));
  
  // Individual costs fit uint64_t while totals, budgets, and positions may exceed it
  const uint64_t maximum=std::numeric_limits<uint64_t>::max();
  const RunOptions wideOptions{.L=1<<20,.T=7,.reservedRngDoublesPerSite=maximum};
  RngLayout wide(wideOptions);
  wide.addOperation(RngScope::Hit,maximum);
  const auto wideOperation=
    wide.addOperation(RngScope::Hit,maximum);
  const auto wideRange=
    wide.operationRange(rng,wideOperation,maximum,maximum-1);
  assert(equalPositions(wide.doubleDrawCosts[1],RngPosition{{{maximum-1,1,0,0}}}));
  assert(equalPositions(wideRange.position,RngPosition{{{8070450532247928831ULL,10376293541461622785ULL,10376293541461622783ULL,8070450532247928831ULL}}}));
  assert(equalPositions(wideRange.rng.state+wideRange.integerOffset(0),
                        RngPosition{{{16140901064495857662ULL,2305843009213693954ULL,2305843009213693951ULL,16140901064495857663ULL}}}));
  assert(equalPositions(wideRange.rng.state+wideRange.integerOffset(maximum-1),
                        RngPosition{{{16140901064495857658ULL,2305843009213693956ULL,2305843009213693951ULL,16140901064495857663ULL}}}));
  assert(wideRange.uniformAt(maximum-1)>=0.0 and wideRange.uniformAt(maximum-1)<1.0);
  
  if(narg>1)
    {
      const std::string_view mode=arg[1];
      
      if(mode=="rng-insufficient-reservation")
        exact.addOperation(RngScope::Hit,1);
      else if(mode=="rng-insufficient-conf-reservation")
        exact.addOperation(RngScope::GaugeConf,1);
      else if(mode=="rng-draw-outside-range")
        origin.uniformAt(NDIM);
      else if(mode=="rng-hit-index-overflow")
        layout.operationRange(rng,layout.origin,3,maximum);
      else if(mode=="rng-empty-operation-draw")
        layout.operationRange(rng,deterministic,3).uniformAt(0);
      else if(mode=="rng-missing-hit-index")
        layout.operationRange(rng,layout.origin,3);
      else if(mode=="rng-unexpected-hit-index")
        layout.operationRange(rng,gauge,3,0);
      else if(mode=="rng-invalid-scope")
        layout.addOperation(static_cast<RngScope>(2),1);
      else if(mode=="rng-allocation-outside-budget")
        exact.operationRange(rng,{RngScope::Hit,RngPosition{{{4,0,0,0}}},1},0,0);
    }
}

/// Verify origin placement and site-fastest operation draws on the configured backend
void validateRngExecution()
{
  using namespace nissa;
  
  RunOptions options{.L=4,.T=4,.nHits=8};
  RngLayout layout(options);
  layout.origin=
    layout.addOperation(RngScope::Hit,options.drawGlobalOriginTask.doubleDrawCost());
  const uint64_t volume=static_cast<uint64_t>(nissa::glbVol);
  const auto sourceOperation=
    layout.addOperation(RngScope::Hit,checkedMultiply(volume,uint64_t{2}));
  const Sitmo::Rng rng(3472291050U);
  const auto originRange=
    layout.operationRange(rng,layout.origin,3,5);
  const auto origin=
    options.drawGlobalOriginTask.execute(originRange);
  
  // The origin's four draws are consecutive, independent of the lattice volume
  for(int mu=0;mu<NDIM;mu++)
    {
      const RngPosition raw{{{3'072'000'000+uint64_t{2}*static_cast<uint64_t>(mu),1'536'000'000,0,0}}};
      const uint32_t low=rng.generateNth(raw);
      const uint32_t high=rng.generateNth(raw+uint64_t{1});
      constexpr double F=4294967296.0;
      const double expected=std::min((low/F+high)/F,0.999999999999999889);
      assert(origin[mu]==static_cast<int>(expected*nissa::glbSize[mu]));
    }
  
  // The actual source callback must place its point at that same global coordinate
  const GaugeRecipes noGaugeRecipes;
  const auto handlers=
    makeNissaTaskHandlers(options,noGaugeRecipes,rng,layout,3,5);
  const auto checkPoint=[](const TaskHandlers& selectedHandlers,
                           const Coords& selectedOrigin)
  {
    LineValue point=selectedHandlers.source(Source{});
    auto shadow=point.field.shadow<nissa::spincolor,nissa::defaultSpaceTimeLayout,
                                    nissa::defaultMemorySpace>("rng point source",
                                      std::vector<size_t>(point.field.outer.size(),0));
    const auto hostPoint=
      shadow.getSurelyReadableOn<nissa::MemorySpace::CPU>();
    
    for(int64_t site=0;site<nissa::locVol;site++)
      {
        bool atOrigin=true;
        
        for(int mu=0;mu<NDIM;mu++)
          atOrigin=atOrigin and nissa::glbCoordOfLoclx[site][mu]==selectedOrigin[mu];
        
        for(int degree=0;degree<NDIRAC*NCOL*2;degree++)
          assert(hostPoint(site,degree)==(atOrigin and degree==0?1.0:0.0));
      }
  };
  checkPoint(handlers,origin);
  
  // CPU and GPU kernels use the global site, irrespective of component evaluation order
  const auto range=
    layout.operationRange(rng,sourceOperation,3,5);
  nissa::LxField<double[2]> field("rng operation validation",nissa::WITHOUT_HALO);
  PAR(0,nissa::locVol,
      CAPTURE(range,volume,TO_WRITE(field)),
      site,
      {
        const uint64_t globalSite=nissa::glblxOfLoclx[site];
        field(site,1)=range.uniformAt(volume+globalSite);
        field(site,0)=range.uniformAt(globalSite);
      });
  const auto host=
    field.getSurelyReadableOn<nissa::MemorySpace::CPU>();
  
  for(int64_t site=0;site<nissa::locVol;site++)
    for(uint64_t component=0;component<2;component++)
      assert(host(site,component)==range.uniformAt(component*volume+nissa::glblxOfLoclx[site]));
  
  // A disabled implicit operation uses the zero origin and owns no random draws
  options.drawGlobalOriginTask.enabled=false;
  RngLayout disabled(options);
  disabled.origin=
    disabled.addOperation(RngScope::Hit,options.drawGlobalOriginTask.doubleDrawCost());
  const auto zero=
    options.drawGlobalOriginTask.execute(disabled.operationRange(rng,disabled.origin,3,5));
  
  for(int mu=0;mu<NDIM;mu++)
    assert(zero[mu]==0);
  
  // A nonzero fixed origin is shared by different hits and gauge configurations
  options.drawGlobalOriginTask.fixedOrigin={3,2,1,0};
  options.validate();
  assert(disabled.origin.doubleDrawCost==0);
  
  for(const size_t iHit : {size_t{0},size_t{7}})
    {
      const auto fixedRange=disabled.operationRange(rng,disabled.origin,iHit,iHit);
      assert(options.drawGlobalOriginTask.execute(fixedRange)==options.drawGlobalOriginTask.fixedOrigin);
      const auto fixedHandlers=
        makeNissaTaskHandlers(options,noGaugeRecipes,rng,disabled,iHit,iHit);
      checkPoint(fixedHandlers,options.drawGlobalOriginTask.fixedOrigin);
    }
}

/// Recover a known physical-basis solution with both Wilson signs and a nonzero clover term
void validatePhysicalPropagator(const bool useExternalSolver=false)
{
  using namespace nissa;
  
  const int savedMultigrid=multiGrid::use_multiGrid;
  multiGrid::use_multiGrid=useExternalSolver?1:0;
  RunOptions options{.L=4,.T=4};
  options.drawGlobalOriginTask.enabled=false;
  RngLayout layout(options);
  layout.origin=layout.addOperation(RngScope::Hit,0);
  GaugeRecipes recipes;
  recipes.push_back({std::make_unique<GaugeField>("propagator test gauge",WITH_HALO_EDGES),{}, {}});
  auto& gauge=*recipes[0].field;
  
  // Coordinate-dependent diagonal SU(3) links produce a genuinely nonzero clover field
  PAR(0,locVol,
      CAPTURE(TO_WRITE(gauge)),
      site,
      {
        quad_su3 value{};
        
        for(int mu=0;mu<NDIM;mu++)
          {
            const double angle=0.13*(mu+1)*(1+glbCoordOfLoclx[site][(mu+1)%NDIM]);
            value[mu][0][0][RE]=cos(angle);
            value[mu][0][0][IM]=sin(angle);
            value[mu][1][1][RE]=cos(angle);
            value[mu][1][1][IM]=-sin(angle);
            value[mu][2][2][RE]=1;
          }
        
        quad_su3_copy(gauge[site],value);
      });
  set_borders_invalid(gauge);
  LxField<spincolor> known("known physical solution",WITH_HALO);
  PAR(0,locVol,
      CAPTURE(TO_WRITE(known)),
      site,
      {
        for(int degree=0;degree<NDIRAC*NCOL*2;degree++)
          known(site,degree)=sin(0.01*(glblxOfLoclx[site]+1)*(degree+1));
      });
  set_borders_invalid(known);
  const double knownNorm=known.norm2();
  const AnisDopPars anis{false,1.0,false};
  const auto handlers=makeNissaTaskHandlers(options,recipes,Sitmo::Rng(7),layout,0,0);
  
  for(const int r : {-1,1})
    for(const double csw : (useExternalSolver?std::array{1.2,0.7}:std::array{0.0,1.2}))
      {
        const Prop prop{.kappa=0.11,.mass=0.17,.r=r,.residue=1e-20,.csw=csw};
        LxField<spincolor> twisted("known twisted solution",WITH_HALO);
        
        // Construct R^dagger=(1-i*r*gamma5)/sqrt(2) directly, independently of Pplus/Pminus
        PAR(0,locVol,
            CAPTURE(r,TO_READ(known),TO_WRITE(twisted)),
            site,
            {
              for(int s=0;s<NDIRAC;s++)
                for(int c=0;c<NCOL;c++)
                  {
                    const double sign=r*base_gamma[5].entr[s][RE];
                    twisted[site][s][c][RE]=(known[site][s][c][RE]+sign*known[site][s][c][IM])/sqrt(2.0);
                    twisted[site][s][c][IM]=(known[site][s][c][IM]-sign*known[site][s][c][RE])/sqrt(2.0);
                  }
            });
        set_borders_invalid(twisted);
        LxField<spincolor> qSource("Q on known solution",WITH_HALO);
        
        if(csw==0.0)
          apply_tmQ(qSource,gauge,prop.kappa,anis,r*prop.mass,twisted);
        else
          {
            LxField<clover_term_t> clover("nonzero test clover");
            clover_term(clover,csw,gauge);
            assert(clover.norm2()>1e-12);
            apply_tmclovQ(qSource,gauge,prop.kappa,anis,clover,r*prop.mass,twisted);
          }
        
        LxField<spincolor> physicalSource("known physical source",WITH_HALO);
        
        // The physical operator is R^dagger*gamma5*Q*R^dagger
        PAR(0,locVol,
            CAPTURE(r,TO_READ(qSource),TO_WRITE(physicalSource)),
            site,
            {
              for(int s=0;s<NDIRAC;s++)
                for(int c=0;c<NCOL;c++)
                  {
                    const double g5=base_gamma[5].entr[s][RE];
                    const double re=g5*qSource[site][s][c][RE];
                    const double im=g5*qSource[site][s][c][IM];
                    physicalSource[site][s][c][RE]=(re+r*g5*im)/sqrt(2.0);
                    physicalSource[site][s][c][IM]=(im-r*g5*re)/sqrt(2.0);
                  }
            });
        set_borders_invalid(physicalSource);
        LineValue input(options);
        auto source=input.field.shadow<spincolor,defaultSpaceTimeLayout,defaultMemorySpace>(
          "physical input block",std::vector<size_t>{0});
        source=physicalSource;
        const std::array<WeightedLineInput,1> inputs{{{1.0,&input}}};
        const LineValue result=handlers.extend(ExtendSelectingTPars{Pars(prop)},inputs);
        auto solution=result.field.shadow<spincolor,defaultSpaceTimeLayout,defaultMemorySpace>(
          "physical output block",std::vector<size_t>{0});
        LxField<spincolor> difference("physical solution error",WITHOUT_HALO);
        difference=solution;
        difference-=known;
        assert(difference.norm2()/knownNorm<1e-12);
        difference=source;
        difference-=physicalSource;
        assert(difference.norm2()==0.0);
      }
  
  multiGrid::use_multiGrid=savedMultigrid;
}

/// Check IB imports across dilution, copy, and hit slots, or the supplied sample fields
void validateILDGLineExecution(const bool sample=false)
{
  using namespace nissa;
  
  std::string directory;
  
  if(sample)
    directory="projects/nice/examples/out";
  else
    {
      if(is_master_rank())
        {
          char temporary[]{"/tmp/nice-ildg-line-XXXXXX"};
          assert(mkdtemp(temporary)!=nullptr);
          directory=temporary;
        }
      
      broadcast(directory);
    }
  
  // Distinct values in each coordinate expose hit, copy, dilution, and layout mistakes
  const auto expected=[](const size_t hit,
                          const size_t copy,
                          const size_t spin,
                          const size_t color,
                          const int64_t site,
                          const int degree)
  {
    return static_cast<double>(1000000*hit+100000*copy+10000*spin+1000*color+
                               24*glblxOfLoclx[site]+degree);
  };
  
  for(int shape=0;shape<(sample?1:4);shape++)
    {
      const bool spin=sample or (shape&1);
      const bool color=sample or (shape&2);
      const int copies=sample?1:(shape%2+1);
      const int hits=sample?1:2;
      const std::string stem=sample?"sourceSOURCE":"sourceTEST";
      const auto program=prepareProgram({{ProgramPart::Kind::Immediate,
        std::format("setGeometry({},{}); setRun(7,{},{}); setDilution({},{}); "
                    "setOrigin(0,0,0,0); src=readILDGLine(\"{}/{}\"); "
                    "corrBlock(\"import\"); corrGammas(\"import\",5,5); "
                    "corrTrace(\"import\",src,src);",
                    glbSize[1],glbSize[0],hits,copies,int(spin),int(color),directory,stem),1}});
      const auto& source=std::get<Source>(program.correlations.instructions.front().operation);
      assert(source.ildgPath and source.tSelect==-1);
      const auto layout=prepareRngLayout(program.options,program.correlations,program.gaugeRecipes);
      assert(layout.origin.doubleDrawCost==0);
      const Sitmo::Rng rng(7);
      
      for(size_t hit=0;hit<static_cast<size_t>(hits);hit++)
        {
          // Write portable SciDAC fixtures through the same library interface used by IB
          if(not sample)
            for(size_t copy=0;copy<static_cast<size_t>(copies);copy++)
              for(size_t s=0;s<(spin?NDIRAC:1);s++)
                for(size_t c=0;c<(color?NCOL:1);c++)
                  {
                    LxField<spincolor> field("ILDG test input",WITHOUT_HALO);
                    PAR(0,locVol,
                        CAPTURE(hit,copy,s,c,TO_WRITE(field)),
                        site,
                        {
                          for(int degree=0;degree<24;degree++)
                            field(site,degree)=1000000*hit+100000*copy+10000*s+1000*c+
                              24*glblxOfLoclx[site]+degree;
                        });
                    const std::string suffix=copies>1?"_copy"+std::to_string(copy):"";
                    const std::string path=directory+"/hit"+std::to_string(hit)+"_"+stem+suffix+
                      "_idso"+std::to_string(s)+"_icso"+std::to_string(c);
                    write_real_vector(path,field,"scidac-binary-data");
                  }
          
          const auto handlers=makeNissaTaskHandlers(program.options,program.gaugeRecipes,rng,layout,0,hit);
          const LineValue line=handlers.source(source);
          assert(line.nCopies==static_cast<size_t>(copies));
          assert(line.diluteSpin==spin and line.diluteColor==color);
          
          for(size_t copy=0;copy<line.nCopies;copy++)
            for(size_t s=0;s<(spin?NDIRAC:1);s++)
              for(size_t c=0;c<(color?NCOL:1);c++)
                {
                  std::vector<size_t> coordinates;
                  
                  if(spin)
                    coordinates.push_back(s);
                  
                  coordinates.push_back(copy);
                  
                  if(color)
                    coordinates.push_back(c);
                  
                  auto shadow=line.field.shadow<spincolor,defaultSpaceTimeLayout,defaultMemorySpace>(
                    "ILDG test output",coordinates);
                  const auto host=shadow.getSurelyReadableOn<MemorySpace::CPU>();
                  double sum=0;
                  double norm=0;
                  
                  for(int64_t site=0;site<locVol;site++)
                    for(int degree=0;degree<24;degree++)
                      {
                        const double value=host(site,degree);
                        
                        if(sample)
                          {
                            if(glbCoordOfLoclx[site][0]==7 and degree==int(2*(s*NCOL+c)))
                              assert(value==1.0 or value==-1.0);
                            else
                              assert(value==0.0);
                            
                            sum+=value;
                            norm+=value*value;
                          }
                        else
                          assert(value==expected(hit,copy,s,c,site,degree));
                      }
                  
                  if(sample)
                    {
                      non_loc_reduce(&sum);
                      non_loc_reduce(&norm);
                      assert(sum==2.0 and norm==64.0);
                    }
                }
        }
    }
  
  if(not sample and is_master_rank())
    {
      std::error_code error;
      std::filesystem::remove_all(directory,error);
      assert(not error);
    }
}

/// Check that a parsed script drives the scheduled source, extensions, and trace for each hit
void validatePreparedExecution()
{
  PreparedProgram program=prepareProgram({{ProgramPart::Kind::Immediate,
    "setGeometry(4,4); setRun(7,2,1); rawRecipe=loadedConf(); corrBlock(\"pion\"); "
    "corrGammas(\"pion\",5,5); source=Source(0); "
    "backward=Prop(.gaugeRecipe=rawRecipe,.kappa=0.12,.r=1)*source; "
    "forward=Gamma(1)*backward; corrTrace(\"pion\",backward,forward);",1}});
  program.configurationEntries.push_back({"","unused","script",1});
  assert(program.gaugeRecipes.size()==1);
  program.gaugeRecipes[0].produce=[](const ConfigurationEntry&,
                                     const GaugeRecipes&,
                                     const OperationRngRange& range)
  {
    assert(range.doubleDrawCost==17);
    auto field=std::make_unique<GaugeField>("script test conf",nissa::WITH_HALO_EDGES);
    nissa::generate_cold_lx_conf(*field);
    
    return field;
  };
  size_t requirementCalls=0;
  program.gaugeRecipes[0].doubleDrawCost=[&requirementCalls](const RunOptions&)
  {
    requirementCalls++;
    
    return uint64_t{17};
  };
  const RngLayout rngLayout=
    prepareRngLayout(program.options,program.correlations,program.gaugeRecipes);
  assert(rngLayout.origin.scope==RngScope::Hit and rngLayout.origin.doubleDrawCost==NDIM);
  assert(rngLayout.gaugeOperations[0].scope==RngScope::GaugeConf);
  assert(rngLayout.gaugeOperations[0].doubleDrawCost==17);
  
  // Disabling origins removes their draws without shifting gauge operations
  RunOptions fixedOrigin=program.options;
  fixedOrigin.drawGlobalOriginTask.enabled=false;
  const RngLayout fixedOriginLayout=
    prepareRngLayout(fixedOrigin,program.correlations,program.gaugeRecipes);
  assert(fixedOriginLayout.origin.doubleDrawCost==0);
  
  // More requested hits preserve earlier positions without an explicit hit-count reservation
  RunOptions moreHits=program.options;
  moreHits.nHits=4;
  const RngLayout moreHitLayout=
    prepareRngLayout(moreHits,program.correlations,program.gaugeRecipes);
  const Sitmo::Rng testRng(7);
  const auto firstOrigin=rngLayout.operationRange(testRng,rngLayout.origin,1,0);
  const auto firstOriginWithMoreHits=moreHitLayout.operationRange(testRng,moreHitLayout.origin,1,0);
  
  for(size_t word=0;word<4;word++)
    assert(firstOrigin.position[word]==firstOriginWithMoreHits.position[word]);
  
  size_t sources=0;
  size_t extensions=0;
  size_t traces=0;
  const TaskHandlerFactory makeHandlers=[&](const RunOptions& options,
                                             const GaugeRecipes& recipes,
                                             const size_t iHit)
  {
    assert(recipes.size()==1 and options.nHits==2 and iHit<2);
    assert(&gaugeField(recipes,0)==recipes[0].field.get());
    
    return TaskHandlers{
      .source=[&,iHit](const Source&)
      {
        sources++;
        LineValue result(options,nissa::MemorySpace::CPU,nissa::SpaceTimeLayout::CPU);
        result.field.tensor.data[0]=static_cast<double>(iHit+1);
        
        return result;
      },
      .extend=[&](const ExtendSelectingTPars& operation,
                  const std::span<const WeightedLineInput> inputs)
      {
        extensions++;
        assert(inputs.size()==1);
        if(const auto* prop=std::get_if<Prop>(&operation.pars))
          assert(prop->gaugeRecipe==0);
        LineValue result(options,nissa::MemorySpace::CPU,nissa::SpaceTimeLayout::CPU);
        result.field.tensor.data[0]=inputs[0].weight*inputs[0].value->field.tensor.data[0]+
          (std::holds_alternative<Prop>(operation.pars)?3.0:5.0);
        
        return result;
      },
      .trace=[&](const TraceContraction& operation,
                 const LineValue& backward,
                 const LineValue& forward)
      {
        traces++;
        const std::set<std::pair<int,int>> expected{{5,5}};
        assert(operation.gammaList==expected);
        Tensor result{{ScalarKind::Complex,{}},nissa::MemorySpace::CPU};
        result.data[0]=backward.field.tensor.data[0]+forward.field.tensor.data[0];
        
        return result;
      }};
  };
  const Sitmo::Rng rng(static_cast<uint32_t>(program.options.seed));
  const ConfigurationResults results=executePreparedConfiguration(program,0,makeHandlers,rng,rngLayout,0);
  assert(requirementCalls==3);
  assert(results.hits.size()==2 and results.entry.sourceList=="script");
  assert(sources==2 and extensions==4 and traces==2);
  assert(program.correlations.outputs.size()==1);
  const auto& output=program.correlations.outputs.front();
  assert(output.name=="pion" and output.traces.size()==1);
  for(const TaskValues& hit : results.hits)
    for(size_t id=0;id<hit.size();id++)
      assert(hit[id].has_value()==(id==output.traces[0]));
  assert(std::get<Tensor>(*results.hits[0][output.traces[0]]).data[0]==13.0);
  assert(std::get<Tensor>(*results.hits[1][output.traces[0]]).data[0]==15.0);
  
  // Gauge-only execution must not invoke the unused per-hit origin operation
  program.correlations=CompiledCorrelations{{},{}};
  program.validate();
  const auto gaugeOnlyLayout=prepareRngLayout(program.options,program.correlations,program.gaugeRecipes);
  const auto gaugeOnly=executePreparedConfiguration(program,0,rng,gaugeOnlyLayout,0);
  assert(gaugeOnly.hits.empty());
}

/// Exercise loading, chained spatial APE, and field release between configuration entries
void validateGaugeExecution()
{
  PreparedProgram prepared=prepareProgram({{ProgramPart::Kind::Immediate,
    "setGeometry(4,4); raw=loadedConf(); zero=apeSmearConf(raw,0.4,0); "
    "smooth=apeSmearConf(raw,0.4,2); chained=apeSmearConf(smooth,0.2,1);",1}});
  nissa::initGrid(prepared.options.T,prepared.options.L);
  assert(nissa::memoryManager(nissa::MemorySpace::CPU)==
         nissa::memoryManager<nissa::MemorySpace::CPU>());
  validateRngExecution();
  validateTaskExecution(prepared.options);
  validatePreparedExecution();
  
  // A line owns one tensor; selecting outer indices borrows a contiguous field block
  TaskValues values(6);
  RunOptions lineOptions=prepared.options;
  lineOptions.nCopies=2;
  lineOptions.diluteSpin=true;
  lineOptions.diluteColor=true;
  values[0]=LineValue{lineOptions,nissa::MemorySpace::CPU,nissa::SpaceTimeLayout::CPU};
  auto& lineValue=std::get<LineValue>(*values[0]);
  auto& line=lineValue.field;
  assert(lineValue.nCopies==2 and lineValue.diluteSpin and lineValue.diluteColor);
  assert(line.outer.size()==3 and line.inner.size()==2);
  assert(line.tensor.descriptor.comps.size()==6);
  const auto selected=line.view(std::array<size_t,3>{1,1,2});
  assert(selected.data==line.tensor.data+11*line.nSites*line.nInner*2);
  const auto typedSelected=line.view({{StandardComp::Color,CompSide::Column,2},
                                      {StandardComp::Spin,CompSide::Column,1},
                                      {StandardComp::Copy,CompSide::None,1}});
  assert(typedSelected.data==selected.data);
  assert(selected.getPtrTo(1,0)==selected.data+line.nInner*2);
  assert(selected.getPtrTo(1,1)==selected.data+(line.nInner+1)*2);
  {
    auto shadow=line.shadow<nissa::spincolor,nissa::SpaceTimeLayout::CPU,
                            nissa::MemorySpace::CPU>("line shadow",std::array<size_t,3>{1,1,2});
    assert(not shadow.ownsData and shadow._data==selected.data);
    auto reference=shadow.getWritable();
    auto shallowCopy=reference;
    assert(shadow.nRef==2);
    assert(shallowCopy.template getPtr<nissa::MemorySpace::CPU>()==selected.data);
  }
  
  // A Lorentz-vector field uses the same owner and exposes one copy at a time
  values[1]=nissa::nice::SField{ScalarKind::Complex,
                                {{StandardComp::Lorentz,CompSide::None,NDIM}},
                                {{StandardComp::Copy,CompSide::None,2}},
                                nissa::MemorySpace::CPU,nissa::SpaceTimeLayout::CPU};
  auto& vectorField=std::get<nissa::nice::SField>(*values[1]);
  const auto vectorView=vectorField.view(std::array<size_t,1>{1});
  assert(vectorView.data==vectorField.tensor.data+vectorField.nSites*NDIM*2);
  vectorView.getPtrTo(0,0)[0]=2.5;
  assert(vectorView.getPtrTo(0,0)[0]==2.5);

  // Labels distinguish repeated outer components of the same kind
  nissa::nice::SField labelled{ScalarKind::Real,{},
                              {{StandardComp::Copy,CompSide::None,2,0},
                               {StandardComp::Copy,CompSide::None,3,1}},
                              nissa::MemorySpace::CPU,nissa::SpaceTimeLayout::CPU};
  const auto labelledView=labelled.view({{StandardComp::Copy,CompSide::None,2,1},
                                        {StandardComp::Copy,CompSide::None,1,0}});
  assert(labelledView.data==labelled.tensor.data+5*labelled.nSites);
  
  // Intern names once; tensor descriptors then compare only numeric kinds
  CompRegistry registry;
  const CompKind customKind=registry.intern("Flavor");
  static_assert(StandardComp::compSize<StandardComp::Spin> ==NDIRAC);
  static_assert(StandardComp::compSize<StandardComp::Color> ==NCOL);
  static_assert(StandardComp::compSize<StandardComp::Lorentz> ==NDIM);
  static_assert(StandardComp::compSize<StandardComp::Site> ==StandardComp::ContextDependentSize);
  static_assert(StandardComp::compRequiresOrientation<StandardComp::Spin>);
  static_assert(StandardComp::compRequiresOrientation<StandardComp::Color>);
  static_assert(not StandardComp::compRequiresOrientation<StandardComp::Site>);
  static_assert(not StandardComp::compRequiresOrientation<StandardComp::Lorentz>);
  static_assert(not StandardComp::compRequiresOrientation<StandardComp::Copy>);
  assert(customKind>StandardComp::Copy);
  assert(registry.intern("Flavor")==customKind);
  assert(registry.intern("Spin")==StandardComp::Spin);
  assert(registry.name(customKind)=="Flavor");
  Comp{customKind,CompSide::Row,2}.validate(&registry);
  values[2]=Tensor{{ScalarKind::Real,
                    {{customKind,CompSide::Row,2},
                     {customKind,CompSide::Row,2,1},
                     {customKind,CompSide::Column,2}}},nissa::MemorySpace::CPU};
  assert(std::get<Tensor>(*values[2]).nDoubles==8);
  
  // An undiluted line has no source spin or color column component
  RunOptions undilutedOptions=prepared.options;
  values[3]=LineValue{undilutedOptions,nissa::MemorySpace::CPU,nissa::SpaceTimeLayout::CPU};
  const auto& undilutedValue=std::get<LineValue>(*values[3]);
  const auto& undiluted=undilutedValue.field;
  assert(undiluted.outer.size()==1 and undiluted.tensor.descriptor.comps.size()==4);
  assert(not lineValue.sameShapeAs(undilutedValue));
  
  // The GPU site layout also works with CPU allocation, so its strides can be checked
  values[4]=LineValue{undilutedOptions,nissa::MemorySpace::CPU,nissa::SpaceTimeLayout::GPU};
  auto& gpuValue=std::get<LineValue>(*values[4]);
  auto& gpuLayout=gpuValue.field;
  assert(gpuValue.sameShapeAs(undilutedValue));
  const auto gpuView=gpuLayout.view(std::array<size_t,1>{0});
  assert(gpuView.getPtrTo(1,0)==gpuView.data+2);
  assert(gpuView.getPtrTo(1,1)==gpuView.data+(gpuLayout.nSites+1)*2);
  
  // An empty gamma selection has a valid shape and no allocated data
  const CompKind gammaKind=registry.intern("GammaPair");
  const CompKind timeKind=registry.intern("TimeSeparation");
  values[5]=Tensor{{ScalarKind::Complex,
                    {{gammaKind,CompSide::None,0},
                     {StandardComp::Copy,CompSide::None,2},
                     {timeKind,CompSide::None,static_cast<size_t>(prepared.options.T)}}},
                   nissa::MemorySpace::CPU};
  const auto& emptyResult=std::get<Tensor>(*values[5]);
  assert(emptyResult.nDoubles==0 and emptyResult.data==nullptr);
  values[0].reset();
  values[1].reset();
  values[2].reset();
  values[3].reset();
  values[4].reset();
  values[5].reset();
  char temporary[]{"/tmp/nice-gauge-XXXXXX"};
  std::string path;
  
  if(nissa::rank==nissa::master_rank)
    {
      const int descriptor=mkstemp(temporary);
      assert(descriptor>=0);
      assert(close(descriptor)==0);
      path=temporary;
    }
  
  nissa::broadcast(path);
  
  // A cold gauge field has unit plaquette and is a fixed point of spatial APE
  {
    GaugeField cold("nice validation cold conf",nissa::WITH_HALO_EDGES);
    nissa::generate_cold_lx_conf(cold);
    nissa::write_ildg_gauge_conf(path,cold);
  }
  
  auto& recipes=prepared.gaugeRecipes;
  int loads=0;
  std::string selectedOutput;
  
  // The first producer must see every slot cleared, including later slots from the previous entry
  auto load=std::move(recipes[0].produce);
  recipes[0].produce=[load=std::move(load),&loads,&selectedOutput](const ConfigurationEntry& entry,
                                                               const GaugeRecipes& current,
                                                               const OperationRngRange& range)
  {
    for(const auto& recipe : current)
      assert(not recipe.field);
    
    loads++;
    selectedOutput=entry.outputDirectory;
    
    return load(entry,current,range);
  };
  
  // The same producers fill their storage again for each configuration entry
  for(const std::string output : {"output first","output second"})
    {
      const RngLayout rngLayout=
        prepareRngLayout(prepared.options,prepared.correlations,recipes);
      prepareGaugeFields(recipes,{path,output},Sitmo::Rng(1),
                         rngLayout,0);
      assert(selectedOutput==output);
      
      for(size_t slot=0;slot<recipes.size();slot++)
        assert(std::abs(nissa::global_plaquette_lx_conf(gaugeField(recipes,slot))-1.0)<1e-12);
    }
  
  assert(loads==2);
  clearGaugeFields(recipes);

  // A lost claim after loading the first field must skip every dependent recipe
  const RngLayout rngLayout=
    prepareRngLayout(prepared.options,prepared.correlations,recipes);
  bool checkedClaim=false;
  const bool claimed=
    prepareGaugeFields(recipes,{path,"lost claim"},Sitmo::Rng(1),
                       rngLayout,0,
                       [&]()
                       {
                         checkedClaim=true;
                         assert(recipes[0].field);
                         for(size_t slot=1;slot<recipes.size();slot++)
                           assert(not recipes[slot].field);

                         return false;
                       });
  assert(not claimed and checkedClaim and loads==3);
  
  for(const auto& recipe : recipes)
    assert(not recipe.field and recipe.produce);
  
  if(nissa::rank==nissa::master_rank)
    assert(std::remove(path.c_str())==0);
}

/// Run frontend checks normally, or a separate initialized-grid gauge execution mode
int main(int narg,
         char** arg)
{
#ifdef USE_DDALPHAAMG
  if(narg>1 and std::string_view(arg[1])=="dd-propagator-execution")
    {
      nissa::initNissa(narg,arg);
      nissa::initGrid(4,4);
      
      // Exercise DDalphaAMG export and reinitialization against a known physical solution
      validatePhysicalPropagator(true);
      DD::finalize();
      nissa::closeNissa();
      
      return 0;
    }
#endif
  
  if(narg>1 and std::string_view(arg[1])=="script-execution")
    {
      nissa::initNissa(narg,arg);
      nissa::initGrid(4,4);
      validateRngExecution();
      validatePhysicalPropagator();
      validateILDGLineExecution();
      validatePreparedExecution();
      nissa::closeNissa();
      
      return 0;
    }
  
  if(narg>1 and std::string_view(arg[1])=="ildg-line-sample")
    {
      nissa::initNissa(narg,arg);
      nissa::initGrid(8,4);
      validateILDGLineExecution(true);
      nissa::closeNissa();
      
      return 0;
    }
  
  if(narg>1 and std::string_view(arg[1])=="gauge-execution")
    {
      nissa::initNissa(narg,arg);
      validateGaugeExecution();
      nissa::closeNissa();
      
      return 0;
    }
  
  nissa::init_MPI_thread(narg,arg);
  nissa::get_MPI_nranks();
  nissa::get_MPI_rank();
  
  // Completing an empty exchange is valid, including when a single rank has no halo traffic
  nissa::waitAsyncCommsFinish({});
  
  if(narg>1 and std::string_view(arg[1])=="program-pieces")
    {
      validateProgramPieces(narg,arg);
      
#ifdef USE_MPI
      MPI_Finalize();
#endif
      
      return 0;
    }
  
  if(narg>1 and std::string_view(arg[1])=="tensor-invalid-spin-extent")
    TensorDescriptor{ScalarKind::Real,
                     {{StandardComp::Spin,CompSide::Row,NDIRAC+1}}}.validate();
  
  if(narg>1 and std::string_view(arg[1])=="tensor-repeated-label")
    TensorDescriptor{ScalarKind::Real,
                     {{StandardComp::Spin,CompSide::Row,NDIRAC},
                      {StandardComp::Spin,CompSide::Row,NDIRAC}}}.validate();
  
  if(narg>1 and std::string_view(arg[1])=="tensor-missing-orientation")
    TensorDescriptor{ScalarKind::Real,
                     {{StandardComp::Spin,CompSide::None,NDIRAC}}}.validate();
  
  if(narg>1 and std::string_view(arg[1])=="tensor-unexpected-orientation")
    TensorDescriptor{ScalarKind::Real,
                     {{StandardComp::Lorentz,CompSide::Row,NDIM}}}.validate();
  
  if(narg>1 and std::string_view(arg[1])=="tensor-named-component")
    {
      CompRegistry registry;
      const CompKind kind=registry.intern("Flavor");
      Comp{kind,static_cast<CompSide>(-1),2}.validate(&registry);
    }
  
  validateRngLayout(narg,arg);
  validateLaunchArgs(narg,arg);
  validateProgramPieces(narg,arg);
  validateConfigurationSources(narg,arg);
  validateRunSettings(narg,arg);
  validateGaugeRecipeBindings(narg,arg);
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
