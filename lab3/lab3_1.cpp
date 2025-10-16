#include <iostream>
#include <string>
using namespace std;

/*
secventa de adn

TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA

*/


double calculateTmBasic(const string& dnaSequence) {
    int g_count = 0;
    int c_count = 0;
    int a_count = 0;
    int t_count = 0;

    for (char nucleotide : dnaSequence) {
        if (nucleotide == 'G' || nucleotide == 'g') {
            g_count++;
        } else if (nucleotide == 'C' || nucleotide == 'c') {
            c_count++;
        } else if (nucleotide == 'A' || nucleotide == 'a') {
            a_count++;
        } else if (nucleotide == 'T' || nucleotide == 't') {
            t_count++;
        }
    }

    //formula: Tm = 4(G + C) + 2(A + T)
    return static_cast<double>(4 * (g_count + c_count) + 2 * (a_count + t_count));
}

int main() {
    string dna_sequence;

    cout << "Enter the DNA sequence: ";
    cin >> dna_sequence;

    double melting_temp = calculateTmBasic(dna_sequence);

    cout << "\n--- Calculation Result (Basic Formula) ---" << endl;
    cout << "For sequence: " << dna_sequence << endl;
    cout << "The calculated Melting Temperature (Tm) is: " << melting_temp << (char)248 << "C" << endl;
    return 0;
}
