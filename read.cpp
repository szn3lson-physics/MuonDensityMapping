#include <iostream>
#include <fstream>
#include <string>
#include <conio.h>

using namespace std;

string read_line(int position, int size, string line);

int main (){
    fstream dane;
    dane.open("C:\\Users\\Kacperakis\\Downloads\\dane_1.log", ios::in);

    if(dane.good()){
        cout << "Uzyskano dostep do pliku!" << endl;
        string linia;

        fstream plik("C:\\Users\\Kacperakis\\Downloads\\angles.txt", ios::out);

        while (!dane.eof()) {
            getline(dane, linia);
            if (read_line(2, linia.size(), linia) == "000000"){
            plik << read_line(12, linia.size(), linia) << endl;
            }
        }
    } 
    else {
        cout << "Dostep do pliku zostal zabroniony!" << endl;
    }

    dane.close();
    return 0;
}

string read_line(int position, int line_size, string line){
    string cell = "";
    int count = 0;
    for (int i = 0; i < line_size; i++) {
        char c = line[i];
        if (c == ';') count++;
        else if (count == position) cell += c;
        else if (count > position) break;
    }
    return cell;
} 