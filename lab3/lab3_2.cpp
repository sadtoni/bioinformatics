/*
use AI to design an application that uses the sliding window methodology
in order to scan a sequence of DNA from a FASTA file, and display the melting temperature
along the sequence, by using a chart. the chart should have 2 signals, one for each formula.
the sliding window should have 9 positions.

*/
#include <iostream>
#include <string>
#include <vector>
#include <cmath> // Required for log10

using namespace std;

// Data structure to hold the results for each window position
struct TmPoint {
    int position;       // Starting position of the window (0-indexed)
    double tm_basic;    // Tm from the basic formula
    double tm_advanced; // Tm from the salt-adjusted formula
};

// --- Assume calculateTmBasic() and calculateTmAdvanced() functions ---
//     (as defined in previous answers) exist here.

/**
 * @brief Scans a DNA sequence using a sliding window to calculate Tm.
 *
 * @param fullSequence The complete DNA sequence string.
 * @param windowSize The size of the sliding window (e.g., 9).
 * @param na_concentration The molar concentration of Na+ ions.
 * @return A vector of TmPoint structs, one for each window position.
 */
vector<TmPoint> analyzeSequenceWithSlidingWindow(const string& fullSequence, int windowSize, double na_concentration) {
    vector<TmPoint> results;

    // Check if the sequence is long enough for at least one window
    if (fullSequence.length() < windowSize) {
        cout << "Sequence is shorter than the window size." << endl;
        return results; // Return empty vector
    }

    // Slide the window across the sequence
    // The loop stops when the window can no longer fit
    for (int i = 0; i <= fullSequence.length() - windowSize; ++i) {
        // Extract the 9-base-pair subsequence
        string subSequence = fullSequence.substr(i, windowSize);

        // Calculate Tm for the subsequence using both formulas
        double basic = calculateTmBasic(subSequence);
        double advanced = calculateTmAdvanced(subSequence, na_concentration);

        // Store the results for this position
        results.push_back({i, basic, advanced});
    }

    return results;
}
