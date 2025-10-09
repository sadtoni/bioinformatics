#include <iostream>
#include <string>
#include <set>
#include <vector>
#include <cctype>
#include <iomanip>

//TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA

using namespace std;

void generateCombinations(const string& alphabet, int length, string current, vector<string>& result) {
    if ((int)current.size() == length) {
        result.push_back(current);
        return;
    }

    for (char c : alphabet) {
        generateCombinations(alphabet, length, current + c, result);
    }
}

int countOccurrences(const string& text, const string& pattern) {
    int count = 0;
    int n = text.size();
    int m = pattern.size();

    for (int i = 0; i <= n - m; ++i) {
        bool match = true;
        for (int j = 0; j < m; ++j) {
            if (toupper(text[i + j]) != toupper(pattern[j])) {
                match = false;
                break;
            }
        }
        if (match) count++;
    }

    return count;
}

int main() {
    string input;
    cout << "Enter a string: ";
    getline(cin, input);

    set<char> uniqueLetters;
    for (char c : input) {
        if (isalpha(c)) uniqueLetters.insert(toupper(c));
    }

    string alphabet(uniqueLetters.begin(), uniqueLetters.end());
    cout << "Alphabet of string: " << alphabet << endl;

    vector<string> combos2, combos3;
    generateCombinations(alphabet, 2, "", combos2);
    generateCombinations(alphabet, 3, "", combos3);

    cout << fixed << setprecision(2);

    int totalPairs = max(0, (int)input.size() - 1);
    cout << "\n--- 2-character combinations ---\n";
    for (const string& combo : combos2) {
        int count = countOccurrences(input, combo);
        double percentage = totalPairs > 0 ? (100.0 * count / totalPairs) : 0.0;
        cout << combo << ": " << percentage << "% (" << count << " occurrences)\n";
    }

    int totalTriplets = max(0, (int)input.size() - 2);
    cout << "\n--- 3-character combinations ---\n";
    for (const string& combo : combos3) {
        int count = countOccurrences(input, combo);
        double percentage = totalTriplets > 0 ? (100.0 * count / totalTriplets) : 0.0;
        cout << combo << ": " << percentage << "% (" << count << " occurrences)\n";
    }

    return 0;
}
