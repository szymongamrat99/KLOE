#include <ErrorLogs.h>

/**
 * @brief Find clusters not associated with any track.
 * @param nclu Number of clusters (Fortran: 1-based, C++: 1..nclu)
 * @param ntcl Number of tracks
 * @param asscl Array of associations (size ntcl), asscl[j] = cluster index (Fortran: 1-based)
 * @param NCLMIN Minimal number of clusters required
 * @param logger Error logger (ErrorHandling::ErrorLogs)
 * @param neuclulist Output: vector of cluster indices (1-based) not associated with any track
 * @return true if error (less than NCLMIN clusters found), false otherwise
 */
int find_neuclu(
    int nclu,
    int ntcl,
    const int* asscl,
    int NCLMIN,
    ErrorHandling::ErrorLogs& logger,
    std::vector<int>& neuclulist
) {
    neuclulist.clear();
    int neucluind = 0;
    for (int i = 1; i <= nclu; ++i) { // Fortran: 1-based
        bool neuclu = true;
        for (int j = 1; j <= ntcl; ++j) {
            if (asscl[j - 1] == i) { // C++: 0-based
                neuclu = false;
                break;
            }
        }
        if (neuclu) {
            ++neucluind;
            neuclulist.push_back(i);
        }
    }
    if (neucluind < NCLMIN) {
        auto err = ErrorHandling::ErrorCodes::RANGE;
        logger.getErrLog(err, "find_neuclu: less than NCLMIN clusters found");
        return int(err);
    }
    return 0;
}
