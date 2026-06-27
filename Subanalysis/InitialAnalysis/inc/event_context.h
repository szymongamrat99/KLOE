#ifndef EVENT_CONTEXT_H
#define EVENT_CONTEXT_H
#include <DataAccessWrapper.h>
#include <GeneratedVariables.h>
#include <charged_mom.h>
#include <kloe_class.h>

using namespace KLOE;

class EventContext
{
  public:
    EventContext(DataAccessWrapper &dataAccess, pm00 &obj, ErrorHandling::ErrorLogs &logger, BaseKinematics &baseKin) : fdataAccess(dataAccess), fobj(obj), flogger(logger), fgenVarClassifier(), fBaseKin(baseKin) {}
    ~EventContext() = default;

    void NewEvent();

    int GetNumberOfVerticesWithTwoTracks();
    ErrorHandling::ErrorCodes FilterNeutralClusters(int minClusters, double minClusterEnergy);

    bool CheckChargedVerticesWithTracksAssignedToClusters();
    bool ReconstructChargedParticles();

    void ProcessMonteCarloTruth(bool isMcEnabled, 
                                int mctruthSignal,
                                std::array<UInt_t, 8>& mctruth_num);

    ErrorHandling::ErrorCodes ReconstructKaonIntoChargedPions();
    ErrorHandling::ErrorCodes ReconstructKSIntoChargedParticles();
    ErrorHandling::ErrorCodes ReconstructKLIntoChargedParticles();
    ErrorHandling::ErrorCodes ReconstructKaonClosestToIP();


    // --- Getters for derived variables ---
    void SetCurrentEventAnalysis(KLOE::ChargedVtxRec<>& eventAnalysis) {
        fCurrentEventAnalysis = &eventAnalysis;
    }
    
    // Zwraca obiekt (używane wewnątrz metod kontekstu lub hipotez)
    KLOE::ChargedVtxRec<>& GetEventAnalysis() { return *fCurrentEventAnalysis; }

  private:
    DataAccessWrapper &fdataAccess;
    pm00 &fobj;
    GeneratedVariables fgenVarClassifier;
    ErrorHandling::ErrorLogs &flogger;
    KLOE::ChargedVtxRec<>* fCurrentEventAnalysis = nullptr;
    BaseKinematics &fBaseKin;

    // === Derived variables ===
    std::unordered_map<Int_t, Int_t> fVtxTracksAssignedNum; ///< Map of vertex index to number of tracks assigned

};

#endif // EVENT_CONTEXT_H