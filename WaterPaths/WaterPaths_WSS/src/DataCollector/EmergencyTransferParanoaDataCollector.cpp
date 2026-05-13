//
// EmergencyTransferParanoaDataCollector.cpp
//

#include <iomanip>
#include <sstream>
#include "EmergencyTransferParanoaDataCollector.h"
#include "../Utils/Constants.h"

using namespace Constants;

EmergencyTransferParanoaDataCollector::EmergencyTransferParanoaDataCollector(
        EmergencyTransferParanoa *policy, unsigned long realization)
        : DataCollector(policy->id, "ParanoaTransfer", realization,
                        EMERGENCY_TRANSFER_PARANOA, NON_INITIALIZED),
          transfer_policy(policy) {}

void EmergencyTransferParanoaDataCollector::collect_data() {
    transferred_volumes.push_back(transfer_policy->getTransferredVolume());
}

// ── Compact CSV (one value + trailing comma) ─────────────────────────────────

string EmergencyTransferParanoaDataCollector::printCompactStringHeader() {
    return "ParanoaTransfer,";
}

string EmergencyTransferParanoaDataCollector::printCompactString(int week) {
    stringstream ss;
    ss << transferred_volumes.at((unsigned int) week) << ",";
    return ss.str();
}

// ── Tabular (fixed-width) ────────────────────────────────────────────────────

string EmergencyTransferParanoaDataCollector::printTabularStringHeaderLine1() {
    stringstream ss;
    ss << setw(COLUMN_WIDTH) << "Paranoa";
    return ss.str();
}

string EmergencyTransferParanoaDataCollector::printTabularStringHeaderLine2() {
    stringstream ss;
    ss << setw(COLUMN_WIDTH) << "Transfer";
    return ss.str();
}

string EmergencyTransferParanoaDataCollector::printTabularString(int week) {
    stringstream ss;
    ss << setw(COLUMN_WIDTH) << setprecision(COLUMN_PRECISION)
       << transferred_volumes.at((unsigned int) week);
    return ss.str();
}
