#include <iostream>
using namespace std;

int main() {
    string input;
    cout << "Enter a string: ";
    getline(cin, input);

    int counts[26] = {0};
    int total = 0;

    // Count letters
    for (int i = 0; i < input.size(); i++) {
        char c = input[i];
        if (c >= 'A' && c <= 'Z') {
            c = c - 'A' + 'a';
        }
        if (c >= 'a' && c <= 'z') {
            counts[c - 'a']++;
            total++;
        }
    }

    cout << "\nLetter frequency (as percentage of total letters):\n";
    for (int i = 0; i < 26; i++) {
        if (counts[i] > 0) {
            double percentage = (counts[i] * 100.0) / total;
            cout << char('a' + i) << " : " << percentage << "%" << endl;
        }
    }

    return 0;
}
