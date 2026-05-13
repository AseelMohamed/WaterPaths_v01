//
// Data collector for the EmergencyTransferParanoa drought mitigation policy.
// Records the weekly volume transferred from Paranoa to TortoSM and outputs
// it as a single CSV column ("ParanoaTransfer").
//

#ifndef WATERPATH_EMERGENCYTRANSFERPARANOADATA_H
#define WATERPATH_EMERGENCYTRANSFERPARANOADATA_H

#include <vector>
#include "Base/DataCollector.h"
#include "../DroughtMitigationInstruments/EmergencyTransferParanoa.h"

class EmergencyTransferParanoaDataCollector : public DataCollector {
public:
    explicit EmergencyTransferParanoaDataCollector(EmergencyTransferParanoa *policy,
                                                   unsigned long realization);

    string printTabularString(int week) override;
    string printCompactString(int week) override;
    string printTabularStringHeaderLine1() override;
    string printTabularStringHeaderLine2() override;
    string printCompactStringHeader() override;
    void collect_data() override;

private:
    EmergencyTransferParanoa *transfer_policy;
    vector<double> transferred_volumes;  // one entry per week (hm³)
};

#endif //WATERPATH_EMERGENCYTRANSFERPARANOADATA_H
