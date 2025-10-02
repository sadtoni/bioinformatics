#include <iostream>
#include <string>
#include <set>
#include <cctype>

using namespace std;

int main() {
    string input;
    cout << "Enter a string: ";
    getline(cin, input);

    set<char> uniqueLetters;

    for (char c : input) {
        if (isalpha(c)) {
            c = toupper(c);
            uniqueLetters.insert(c);
        }
    }

    string alphabet;
    for (char c : uniqueLetters) {
        alphabet.push_back(c);
    }

    cout << "Alphabet of string: " << alphabet << endl;
    return 0;
}
