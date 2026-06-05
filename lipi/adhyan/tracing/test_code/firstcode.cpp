#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main()
{
    vector<string> msg {"A","B453","CC"};

    for (const string& word: msg)
    {
        cout << word << " ";
    }
    cout << endl;
}
