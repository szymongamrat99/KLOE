#ifndef EVENT_CONTEXT_H
#define EVENT_CONTEXT_H
#include <DataAccessWrapper.h>
#include <GeneratedVariables.h>
#include <kloe_class.h>

using namespace KLOE;

class EventContext
{
  public:
    EventContext(DataAccessWrapper &dataAccess, pm00 &obj, ErrorHandling::ErrorLogs &logger) : fdataAccess(dataAccess), fobj(obj), flogger(logger), fgenVarClassifier() {}
    ~EventContext() = default;

    void NewEvent();

    int GetNumberOfVerticesWithTwoTracks();
    ErrorHandling::ErrorCodes FilterNeutralClusters(int minClusters, double minEnergy, std::vector<Int_t>& outList);

    bool CheckChargedVerticesWithTracksAssignedToClusters();
    bool ReconstructChargedParticles();

  private:
    DataAccessWrapper &fdataAccess;
    pm00 &fobj;
    GeneratedVariables fgenVarClassifier;
    ErrorHandling::ErrorLogs &flogger;

    std::vector<Int_t> fIv;        ///< Indices of vertices
    std::vector<Int_t> fAssClu;    ///< Indices of clusters assigned to charged tracks
    std::vector<Int_t> fAssTr;     ///< Indices of charged tracks assigned to clusters
    std::vector<Int_t> fNeuCluList; ///< Indices of neutral clusters
    Int_t fNClu;                   ///< Number of clusters

    // === Derived variables ===
    std::unordered_map<Int_t, Int_t> fVtxTracksAssignedNum; ///< Map of vertex index to number of tracks assigned

};

#endif // EVENT_CONTEXT_H