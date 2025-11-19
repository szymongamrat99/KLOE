#include "StatisticalCutter.h"
#include <fstream>
#include <iomanip>
#include <ctime>
#include <sstream>

using json = nlohmann::json;

void StatisticalCutter::GenerateReport(const std::string& reportConfigPath, const std::string& outputDir) const {
    // Wczytaj konfigurację raportu z JSON
    json reportConfig;
    {
        std::ifstream f(reportConfigPath);
        if (!f.is_open()) {
            throw std::runtime_error("Cannot open report config: " + reportConfigPath);
        }
        f >> reportConfig;
    }

    // Pobierz kolumny do raportu
    std::vector<std::string> columns;
    if (reportConfig.contains("reportColumns") && reportConfig["reportColumns"].is_array()) {
        for (const auto& col : reportConfig["reportColumns"]) {
            columns.push_back(col.get<std::string>());
        }
    } else {
        // Domyślne kolumny
        columns = {"cutId", "survived_signal", "survived_background", "efficiency", "purity", "signal_to_bg"};
    }

    // Mapowanie nazwy hipotezy
    std::string hypothesisName;
    switch(hypoCode_) {
        case KLOE::HypothesisCode::SIGNAL:
            hypothesisName = "SIGNAL";
            break;
        case KLOE::HypothesisCode::FOUR_PI:
            hypothesisName = "FOUR_PI";
            break;
        case KLOE::HypothesisCode::OMEGAPI:
            hypothesisName = "OMEGA_PI";
            break;
        case KLOE::HypothesisCode::SIMONA_ANALYSIS:
            hypothesisName = "SIMONA_ANALYSIS";
            break;
        default:
            hypothesisName = "UNKNOWN";
    }

    // Funkcja do generowania jednego raportu
    auto generateSingleReport = [&](bool useFV, const std::string& suffix) {
        std::string filename = outputDir + "/report_" + hypothesisName + suffix + ".csv";
        std::ofstream report(filename);
        
        if (!report.is_open()) {
            throw std::runtime_error("Cannot create report file: " + filename);
        }

        // ==== HEADER - Tekst ====
        report << "# CUTS ANALYSIS REPORT\n";
        
        // Data i czas
        auto now = std::time(nullptr);
        auto tm = *std::localtime(&now);
        report << "# Generated: " << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << "\n";
        report << "#\n";

        // Założenia
        report << "# ASSUMPTIONS:\n";
        report << "#   Hypothesis: " << hypothesisName << "\n";
        report << "#   Fiducial Volume: " << (useFV ? "YES" : "NO") << "\n";
        report << "#   Normalization Mode: " << (normalizationMode_ == NormalizationMode::FIDUCIAL_VOLUME ? "FIDUCIAL_VOLUME" : "TOTAL_EVENTS") << "\n";
        report << "#\n";

        // Statystyki początkowe
        report << "# INITIAL STATISTICS:\n";
        if (useFV) {
            report << "#   Total Signal (in FV): " << totalSignalInFV_ << "\n";
            report << "#   Total Background (in FV): " << totalBackgroundInFV_ << "\n";
            report << "#   Total Events (in FV): " << (totalSignalInFV_ + totalBackgroundInFV_) << "\n";
            report << "#   Total Signal (in FV, excluding errors): " << totalSignalInFVExcludingMinus1_ << "\n";
        } else {
            report << "#   Total Signal: " << totalSignal_ << "\n";
            report << "#   Total Background: " << totalBackground_ << "\n";
            report << "#   Total Events: " << (totalSignal_ + totalBackground_) << "\n";
            report << "#   Total Signal (excluding errors): " << totalSignalExcludingMinus1_ << "\n";
        }
        report << "#\n";

        // Cięcia aktywne - wybierz tylko te które powinny być w raporcie
        // (pomiń członków grup - one liczą się poprzez grupę syntetyczną)
        report << "# ACTIVE CUTS:\n";
        std::vector<size_t> cutsForReport;
        std::set<std::string> seenChannels;  // Śledzenie jakie kanały (grupy) już widzieliśmy
        
        for (size_t idx : activeCutIndices_) {
            if (idx >= cuts_.size()) continue;
            
            const auto& cut = cuts_[idx];
            
            // Jeśli to członek grupy - znajdź syntetyczne cięcie dla jego kanału
            if (groupMemberIndices_.count(idx)) {
                // To członek grupy - sprawdź czy już dodaliśmy syntetyczne cięcie dla tego kanału
                std::string channel = cut.channel;
                if (!seenChannels.count(channel)) {
                    // Znajdź syntetyczne cięcie dla tego kanału
                    for (size_t i = 0; i < cuts_.size(); ++i) {
                        if (cuts_[i].isSyntheticGroup && cuts_[i].channel == channel) {
                            cutsForReport.push_back(i);
                            seenChannels.insert(channel);
                            break;
                        }
                    }
                }
            } else {
                // To normalne cięcie (nie-członek grupy) - dodaj bezpośrednio
                cutsForReport.push_back(idx);
            }
        }
        
        // Wypisz na ekran
        if (cutsForReport.empty()) {
            report << "#   None (all events pass)\n";
        } else {
            for (size_t idx : cutsForReport) {
                if (cuts_[idx].isSyntheticGroup) {
                    report << "#   - " << cuts_[idx].cutId << " (Group)\n";
                } else {
                    report << "#   - " << cuts_[idx].cutId << "\n";
                }
            }
        }
        report << "#\n";

        // Definicje grup (tylko syntetyczne)
        report << "# GROUP DEFINITIONS:\n";
        int groupNum = 1;
        for (size_t idx : cutsForReport) {
            if (cuts_[idx].isSyntheticGroup) {
                report << "# Group" << groupNum << "=" << cuts_[idx].cutId << "\n";
                groupNum++;
            }
        }
        report << "#\n";

        // CSV Header
        report << "CutID";
        for (const auto& col : columns) {
            if (col != "cutId") {
                if (col == "survived_signal") {
                    report << ",Signal";
                } else if (col == "survived_background") {
                    report << ",Background";
                } else if (col == "efficiency") {
                    report << ",Efficiency(%)";
                } else if (col == "efficiency_error") {
                    report << ",Eff.Error(%)";
                } else if (col == "purity") {
                    report << ",Purity(%)";
                } else if (col == "purity_error") {
                    report << ",Pur.Error(%)";
                } else if (col == "signal_to_bg") {
                    report << ",S/B";
                }
            }
        }
        // Dodatkowe kolumny
        report << ",Signal(no-errors),Efficiency-noErr(%),Rel.Eff(%)\n";

        // CSV Data - WSZYSTKIE AKTYWNE CIĘCIA (grupy + normalne)
        double prevSigNoErr = 0.0;  // Do obliczenia wydajności względnej
        bool isFirstCut = true;
        double totalNoErrForFirst = useFV ? totalSignalInFVExcludingMinus1_ : totalSignalExcludingMinus1_;
        
        for (size_t cutIdx : cutsForReport) {
            if (cutIdx >= cuts_.size()) continue;

            const auto& cut = cuts_[cutIdx];
            
            // Liczenie dla tego cięcia
            size_t sig_count = useFV ? survivedSignalInFV_[cutIdx] : survivedSignal_[cutIdx];
            size_t bg_count = useFV ? survivedBackgroundInFV_[cutIdx] : survivedBackground_[cutIdx];
            
            // Liczenie bez błędów (excludingMinus1)
            size_t sig_noErr = useFV ? survivedSignalInFVExcludingMinus1_[cutIdx] : survivedSignalExcludingMinus1_[cutIdx];
            double totalNoErr = useFV ? totalSignalInFVExcludingMinus1_ : totalSignalExcludingMinus1_;

            report << cut.cutId;
            for (const auto& col : columns) {
                if (col != "cutId") {
                    if (col == "survived_signal") {
                        report << "," << sig_count;
                    } else if (col == "survived_background") {
                        report << "," << bg_count;
                    } else if (col == "efficiency") {
                        double total = useFV ? totalSignalInFV_ : totalSignal_;
                        double eff = total > 0 ? (static_cast<double>(sig_count) / total * 100.0) : 0.0;
                        report << "," << std::fixed << std::setprecision(2) << eff;
                    } else if (col == "efficiency_error") {
                        double total = useFV ? totalSignalInFV_ : totalSignal_;
                        double eff = total > 0 ? static_cast<double>(sig_count) / total : 0.0;
                        double effError = total > 0 ? (std::sqrt(eff * (1.0 - eff) / total) * 100.0) : 0.0;
                        report << "," << std::fixed << std::setprecision(2) << effError;
                    } else if (col == "purity") {
                        double purity = (sig_count + bg_count) > 0 ? (static_cast<double>(sig_count) / (sig_count + bg_count) * 100.0) : 0.0;
                        report << "," << std::fixed << std::setprecision(2) << purity;
                    } else if (col == "purity_error") {
                        double purity = (sig_count + bg_count) > 0 ? static_cast<double>(sig_count) / (sig_count + bg_count) : 0.0;
                        double purityError = 0.0;
                        if (sig_count + bg_count > 0 && sig_count > 0) {
                            purityError = purity * std::sqrt(
                                (1.0 - purity) * (1.0 - purity) / sig_count + 
                                purity * purity / (sig_count + bg_count)
                            ) * 100.0;
                        }
                        report << "," << std::fixed << std::setprecision(2) << purityError;
                    } else if (col == "signal_to_bg") {
                        double ratio = bg_count > 0 ? static_cast<double>(sig_count) / bg_count : (sig_count > 0 ? 1e9 : 0.0);
                        report << "," << std::fixed << std::setprecision(3) << ratio;
                    }
                }
            }
            
            // Dodatkowe kolumny: Signal(no-errors), Efficiency-noErr(%), Rel.Eff(%)
            report << "," << sig_noErr;
            
            // Efficiency excluding errors
            double effNoErr = totalNoErr > 0 ? (static_cast<double>(sig_noErr) / totalNoErr * 100.0) : 0.0;
            report << "," << std::fixed << std::setprecision(2) << effNoErr;
            
            // Relative efficiency
            double relEff = 0.0;
            if (isFirstCut) {
                // Pierwsze cięcie: względem totalSignalNoErr (excludingMinus1)
                if (totalNoErrForFirst > 0.0) {
                    relEff = (static_cast<double>(sig_noErr) / totalNoErrForFirst * 100.0);
                }
                isFirstCut = false;
            } else {
                // Następne cięcia: względem poprzedniego cięcia
                if (prevSigNoErr > 0.0) {
                    relEff = (static_cast<double>(sig_noErr) / prevSigNoErr * 100.0);
                }
            }
            report << "," << std::fixed << std::setprecision(2) << relEff;
            
            report << "\n";
            
            // Zapisz dla następnego cięcia
            prevSigNoErr = static_cast<double>(sig_noErr);
        }

        report.close();
    };

    // Generuj oba raporty
    generateSingleReport(false, "_total");
    generateSingleReport(true, "_fiducial_volume");
}
