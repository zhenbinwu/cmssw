#ifndef FWCore_Framework_SubProcess_h
#define FWCore_Framework_SubProcess_h

#include "DataFormats/Provenance/interface/BranchID.h"
#include "DataFormats/Provenance/interface/ProductDescriptionFwd.h"
#include "FWCore/Common/interface/FWCoreCommonFwd.h"
#include "FWCore/Framework/interface/EventSetupProvider.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/PrincipalCache.h"
#include "FWCore/Framework/interface/Schedule.h"
#include "FWCore/Framework/interface/TriggerResultsBasedEventSelector.h"
#include "FWCore/Framework/interface/ProductSelectorRules.h"
#include "FWCore/Framework/interface/ProductSelector.h"
#include "FWCore/ServiceRegistry/interface/ProcessContext.h"
#include "FWCore/ServiceRegistry/interface/ServiceLegacy.h"
#include "FWCore/ServiceRegistry/interface/ServiceRegistry.h"
#include "FWCore/ServiceRegistry/interface/ServiceToken.h"
#include "FWCore/Utilities/interface/Algorithms.h"
#include "FWCore/Utilities/interface/BranchType.h"
#include "FWCore/Utilities/interface/get_underlying_safe.h"
#include "FWCore/Utilities/interface/propagate_const.h"

#include "DataFormats/Provenance/interface/SelectedProducts.h"

#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <vector>

namespace edm {
  class ActivityRegistry;
  class BranchIDListHelper;
  class EventPrincipal;
  class EventSetupImpl;
  class ExceptionCollector;
  class HistoryAppender;
  class LuminosityBlockPrincipal;
  class LumiTransitionInfo;
  class MergeableRunProductMetadata;
  class ModuleTypeResolverMaker;
  class ParameterSet;
  class PathsAndConsumesOfModules;
  class Principal;
  class ProcessBlockTransitionInfo;
  class ProductRegistry;
  class PreallocationConfiguration;
  class RunTransitionInfo;
  class ThinnedAssociationsHelper;
  class SubProcessParentageHelper;
  class WaitingTaskHolder;

  namespace eventsetup {
    class EventSetupsController;
  }
  class SubProcess : public EDConsumerBase {
  public:
    SubProcess(ParameterSet& parameterSet,
               ParameterSet const& topLevelParameterSet,
               std::shared_ptr<ProductRegistry const> parentProductRegistry,
               std::shared_ptr<BranchIDListHelper const> parentBranchIDListHelper,
               ProcessBlockHelperBase const& parentProcessBlockHelper,
               ThinnedAssociationsHelper const& parentThinnedAssociationsHelper,
               SubProcessParentageHelper const& parentSubProcessParentageHelper,
               eventsetup::EventSetupsController& esController,
               ActivityRegistry& parentActReg,
               ServiceToken const& token,
               serviceregistry::ServiceLegacy iLegacy,
               PreallocationConfiguration const& preallocConfig,
               ProcessContext const* parentProcessContext,
               ModuleTypeResolverMaker const* typeResolverMaker);

    ~SubProcess() override;

    SubProcess(SubProcess const&) = delete;             // Disallow copying
    SubProcess& operator=(SubProcess const&) = delete;  // Disallow copying
    SubProcess(SubProcess&&) = default;                 // Allow Moving
    SubProcess& operator=(SubProcess&&) = delete;       // Move not supported by PrincipalCache

    //From OutputModule
    void selectProducts(ProductRegistry const& preg,
                        ThinnedAssociationsHelper const& parentThinnedAssociationsHelper,
                        std::map<BranchID, bool>& keepAssociation);

    SelectedProductsForBranchType const& keptProducts() const { return keptProducts_; }

    // Returns the set of modules whose products may be consumed by
    // modules in this SubProcess or its child SubProcesses
    std::vector<ModuleProcessName> keepOnlyConsumedUnscheduledModules(bool deleteModules);

    void doBeginJob();
    void doEndJob(ExceptionCollector&);

    void doEventAsync(WaitingTaskHolder iHolder,
                      EventPrincipal const& principal,
                      std::vector<std::shared_ptr<const EventSetupImpl>> const*);

    template <typename Traits>
    void doBeginProcessBlockAsync(WaitingTaskHolder iHolder,
                                  ProcessBlockTransitionInfo const& iTransitionInfo,
                                  bool cleaningUpAfterException);

    void doEndProcessBlockAsync(WaitingTaskHolder iHolder,
                                ProcessBlockTransitionInfo const& iTransitionInfo,
                                bool cleaningUpAfterException);

    void doBeginRunAsync(WaitingTaskHolder iHolder, RunTransitionInfo const& iTransitionInfo);

    void doEndRunAsync(WaitingTaskHolder iHolder,
                       RunTransitionInfo const& iTransitionInfo,
                       bool cleaningUpAfterException);

    void doBeginLuminosityBlockAsync(WaitingTaskHolder iHolder, LumiTransitionInfo const& iTransitionInfo);

    void doEndLuminosityBlockAsync(WaitingTaskHolder iHolder,
                                   LumiTransitionInfo const& iTransitionInfo,
                                   bool cleaningUpAfterException);

    void doBeginStream(unsigned int streamID);
    void doEndStream(unsigned int streamID, ExceptionCollector& collector, std::mutex& collectorMutex) noexcept;
    void doStreamBeginRunAsync(WaitingTaskHolder iHolder, unsigned int iID, RunTransitionInfo const&);

    void doStreamEndRunAsync(WaitingTaskHolder iHolder,
                             unsigned int iID,
                             RunTransitionInfo const&,
                             bool cleaningUpAfterException);

    void doStreamBeginLuminosityBlockAsync(WaitingTaskHolder iHolder, unsigned int iID, LumiTransitionInfo const&);

    void doStreamEndLuminosityBlockAsync(WaitingTaskHolder iHolder,
                                         unsigned int iID,
                                         LumiTransitionInfo const&,
                                         bool cleaningUpAfterException);

    void writeLumiAsync(WaitingTaskHolder, LuminosityBlockPrincipal&);

    void clearLumiPrincipal(LuminosityBlockPrincipal&);

    using ProcessBlockType = PrincipalCache::ProcessBlockType;
    void writeProcessBlockAsync(edm::WaitingTaskHolder task, ProcessBlockType);

    void writeRunAsync(WaitingTaskHolder, RunPrincipal const&, MergeableRunProductMetadata const*);

    void clearRunPrincipal(RunPrincipal&);

    void clearProcessBlockPrincipal(ProcessBlockType);

    // Call closeFile() on all OutputModules.
    void closeOutputFiles() {
      ServiceRegistry::Operate operate(serviceToken_);
      schedule_->closeOutputFiles();
      for_all(subProcesses_, [](auto& subProcess) { subProcess.closeOutputFiles(); });
    }

    // Call openFiles() on all OutputModules
    void openOutputFiles(FileBlock& fb) {
      ServiceRegistry::Operate operate(serviceToken_);
      schedule_->openOutputFiles(fb);
      for_all(subProcesses_, [&fb](auto& subProcess) { subProcess.openOutputFiles(fb); });
    }

    void updateBranchIDListHelper(BranchIDLists const&);

    // Call respondToOpenInputFile() on all Modules
    void respondToOpenInputFile(FileBlock const& fb);

    // Call respondToCloseInputFile() on all Modules
    void respondToCloseInputFile(FileBlock const& fb) {
      ServiceRegistry::Operate operate(serviceToken_);
      schedule_->respondToCloseInputFile(fb);
      for_all(subProcesses_, [&fb](auto& subProcess) { subProcess.respondToCloseInputFile(fb); });
    }

    // Call shouldWeCloseFile() on all OutputModules.
    bool shouldWeCloseOutput() const {
      ServiceRegistry::Operate operate(serviceToken_);
      if (schedule_->shouldWeCloseOutput()) {
        return true;
      }
      for (auto const& subProcess : subProcesses_) {
        if (subProcess.shouldWeCloseOutput()) {
          return true;
        }
      }
      return false;
    }

    /// Return a vector allowing const access to all the ModuleDescriptions for this SubProcess

    /// *** N.B. *** Ownership of the ModuleDescriptions is *not*
    /// *** passed to the caller. Do not call delete on these
    /// *** pointers!
    std::vector<ModuleDescription const*> getAllModuleDescriptions() const;

    /// Return the number of events this SubProcess has tried to process
    /// (inclues both successes and failures, including failures due
    /// to exceptions during processing).
    int totalEvents() const { return schedule_->totalEvents(); }

    /// Return the number of events which have been passed by one or more trigger paths.
    int totalEventsPassed() const {
      ServiceRegistry::Operate operate(serviceToken_);
      return schedule_->totalEventsPassed();
    }

    /// Return the number of events that have not passed any trigger.
    /// (N.B. totalEventsFailed() + totalEventsPassed() == totalEvents()
    int totalEventsFailed() const {
      ServiceRegistry::Operate operate(serviceToken_);
      return schedule_->totalEventsFailed();
    }

    /// Return the trigger report information on paths,
    /// modules-in-path, modules-in-endpath, and modules.
    void getTriggerReport(TriggerReport& rep) const {
      ServiceRegistry::Operate operate(serviceToken_);
      schedule_->getTriggerReport(rep);
    }

    /// Return whether each output module has reached its maximum count.
    /// If there is a subprocess, get this information from the subprocess.
    bool terminate() const {
      ServiceRegistry::Operate operate(serviceToken_);
      if (schedule_->terminate()) {
        return true;
      }
      for (auto const& subProcess : subProcesses_) {
        if (subProcess.terminate()) {
          return true;
        }
      }
      return false;
    }

    ///  Clear all the counters in the trigger report.
    void clearCounters() {
      ServiceRegistry::Operate operate(serviceToken_);
      schedule_->clearCounters();
      for_all(subProcesses_, [](auto& subProcess) { subProcess.clearCounters(); });
    }

    void initializePathsAndConsumes();

  private:
    void beginJob();
    void endJob(ExceptionCollector&);
    void processAsync(WaitingTaskHolder iHolder,
                      EventPrincipal const& e,
                      std::vector<std::shared_ptr<const EventSetupImpl>> const*);

    void propagateProducts(BranchType type, Principal const& parentPrincipal, Principal& principal) const;
    bool parentProducedProductIsKept(Principal const& parentPrincipal, Principal& principal) const;
    void fixBranchIDListsForEDAliases(
        std::map<BranchID::value_type, BranchID::value_type> const& droppedBranchIDToKeptBranchID);
    void keepThisBranch(ProductDescription const& desc,
                        std::map<BranchID, ProductDescription const*>& trueBranchIDToKeptBranchDesc,
                        std::set<BranchID>& keptProductsInEvent);

    std::map<BranchID::value_type, BranchID::value_type> const& droppedBranchIDToKeptBranchID() {
      return droppedBranchIDToKeptBranchID_;
    }

    std::shared_ptr<BranchIDListHelper const> branchIDListHelper() const {
      return get_underlying_safe(branchIDListHelper_);
    }
    std::shared_ptr<BranchIDListHelper>& branchIDListHelper() { return get_underlying_safe(branchIDListHelper_); }
    std::shared_ptr<ThinnedAssociationsHelper const> thinnedAssociationsHelper() const {
      return get_underlying_safe(thinnedAssociationsHelper_);
    }
    std::shared_ptr<ThinnedAssociationsHelper> thinnedAssociationsHelper() {
      return get_underlying_safe(thinnedAssociationsHelper_);
    }

    std::shared_ptr<ActivityRegistry> actReg_;  // We do not use propagate_const because the registry itself is mutable.
    ServiceToken serviceToken_;
    std::shared_ptr<ProductRegistry const> parentPreg_;
    std::shared_ptr<ProductRegistry const> preg_;
    edm::propagate_const<std::shared_ptr<BranchIDListHelper>> branchIDListHelper_;
    edm::propagate_const<std::shared_ptr<SubProcessBlockHelper>> processBlockHelper_;
    edm::propagate_const<std::shared_ptr<ThinnedAssociationsHelper>> thinnedAssociationsHelper_;
    edm::propagate_const<std::shared_ptr<SubProcessParentageHelper>> subProcessParentageHelper_;
    std::unique_ptr<ExceptionToActionTable const> act_table_;
    std::shared_ptr<ProcessConfiguration const> processConfiguration_;
    ProcessContext processContext_;
    std::unique_ptr<PathsAndConsumesOfModules> pathsAndConsumesOfModules_;
    //We require 1 history for each Run, Lumi and Stream
    // The vectors first hold Stream info, then Lumi then Run
    unsigned int historyLumiOffset_;
    unsigned int historyRunOffset_;
    std::vector<ProcessHistoryRegistry> processHistoryRegistries_;
    std::vector<HistoryAppender> historyAppenders_;
    PrincipalCache principalCache_;
    //vector index is principal's index value
    std::vector<std::shared_ptr<RunPrincipal>> inUseRunPrincipals_;
    std::vector<std::shared_ptr<LuminosityBlockPrincipal>> inUseLumiPrincipals_;
    edm::propagate_const<std::shared_ptr<eventsetup::EventSetupProvider>> esp_;
    edm::propagate_const<std::unique_ptr<Schedule>> schedule_;
    std::vector<SubProcess> subProcesses_;
    edm::propagate_const<std::unique_ptr<ParameterSet>> processParameterSet_;

    // keptProducts_ are pointers to the ProductDescription objects describing
    // the branches we are to write.
    //
    // We do not own the ProductDescriptions to which we point.
    SelectedProductsForBranchType keptProducts_;
    ProductSelectorRules productSelectorRules_;
    ProductSelector productSelector_;

    //EventSelection
    bool wantAllEvents_;
    ParameterSetID selector_config_id_;
    mutable detail::TriggerResultsBasedEventSelector selectors_;

    // needed because of possible EDAliases.
    // filled in only if key and value are different.
    std::map<BranchID::value_type, BranchID::value_type> droppedBranchIDToKeptBranchID_;
  };

  // free function
  std::vector<ParameterSet> popSubProcessVParameterSet(ParameterSet& parameterSet);
}  // namespace edm
#endif
