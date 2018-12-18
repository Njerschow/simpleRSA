/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: nikitajerschow
 *
 * Created on November 28, 2018, 1:53 AM
 */

#include <cstdlib>
#include <iostream>
#include <bitset>
#include <tuple>
#include <map>
#include <math.h>

namespace constant {
    const int MIL_RAB_ITER = 20;
    const int MAX_EUCL_ITER = 999; // maximum iterations for euclidian algorithm
}

using namespace std;

unsigned char randomBit(unsigned int i, bool print = false) {
    int x = rand();
    if (print) {
        cout << "b_" << i << "|" << bitset<32>(x) << "|" << (x & 0b1) << endl;
    }
    return (x & 0b1);
}

unsigned char binToInt(char binaryString[]) {
    int tmpInt = 0;
    while (*binaryString) {
        tmpInt <<= 1;
        if (*binaryString == '1')
            tmpInt += 1;
        binaryString++;
    }
    return (unsigned char) tmpInt;
}

unsigned int randomNumber(bool print = false) {
    unsigned char c = 0b01000001;
    for (int i = 5; i > 0; i--) {
        unsigned char b = randomBit(i, print);
        c += b << i;
    }
    return (unsigned int) c;
}

unsigned int modulo(int i, int n) {
    if (n == 0) {
        cerr << "ERROR: modulo cannot be zero" << endl;
    }
    while (i >= n) { //i greater than n
        i -= n;
    }
    while (i < 0) { //i negative

        i += n;
    }
    return i;
}

bool primeTest(unsigned int a, unsigned int n, bool print = false) //a^(n-1) mod n
{
    if (print == true) {
        cout << "n = " << n << ", a = " << a << endl;
        cout << "i\t|xi\t|z\t|y\t|y" << endl;
    }

    int y = 1;

    bitset<32> bs = bitset<32>(n - 1);
    size_t i = bs.size();
    while (bs[i - 1] != 1) { //skip logic to skip padding of the binary number
        i--;
    }

    for (i; i > 0; i--) {
        int x = bs[i - 1];
        // minus one because I started at bs.size(). 
        //if the for loop were i=bs.size() -1; i>=0; i--, the stop condition would never be satisfied due to underflow

        int z = y;
        y = modulo(y*y, n);

        if (print) cout << i << "\t|" << x << "\t|" << z << "\t|" << y << "\t|"; //TODO: check printing in this func

        if (y == 1 && z != 1 && z != n - 1) {
            if (print) cout << n << " is not prime because " << y << "^2 mod " << n << " = 1 and " << y << " != 1 and" << y << " != " << n << " - 1" << endl;
            return false;
        }
        if (x == 1) {
            y = modulo(y*a, n);
        }

        if (print) cout << y << endl;
    }
    if (y != 1) {
        if (print) cout << n << " is not prime because " << a << "^(" << n << "-1) mod " << n << " is not equal to 1" << endl;
        return false;
    } else {
        return true;
    }
}

bool isPrime(unsigned int n, bool print = false) {
    bool isPrime = true;
    for (int i = 0; i < constant::MIL_RAB_ITER; i++) {
        unsigned int a = modulo(rand(), n);
        if (a <= 1) { //if a is zero or one, discard and retry
            i--;
            continue;
        }
        isPrime = primeTest(a, n);

        if (!isPrime) {
            if (print) primeTest(a, n, true); // run again to print
            return isPrime;
        }
    }
    if (isPrime && print) isPrime = primeTest(modulo(rand(), n), n, true); // do printing required, I realize this might be a composite, but with very low probability
    return isPrime;
}

unsigned int getPrime() {
    unsigned int potPrime = randomNumber();
    while (!isPrime(potPrime)) { // ensure we have a prime
        potPrime = randomNumber();
    }
    return potPrime;
}

unsigned int fastExp(unsigned int a, unsigned int x, unsigned int n, bool print = false) {
    unsigned int y = 1;
    bitset<32> bx = bitset<32>(x);
    size_t i = bx.size();

    if (print) cout << "i\t|xi\t|y\t|y" << endl;

    while (bx[i - 1] != 1) { //skip logic to skip padding of the binary number
        i--;
    }
    for (; i > 0; i--) {
        int xi = bx[i - 1]; // because of indexing issues as described earlier

        y = modulo(y*y, n);

        if (print) cout << i << "\t|" << xi << "\t|" << y << "\t|";
        if (xi == 1) {
            y = modulo(y*a, n);
        }

        if (print) cout << y << endl;
    }
    return y;
}

/* The ExtEuclidian function performs the extended euclidian algorithm to find the gcd of e and phi(n)
 * and the multiplicative inverse of e mod phi(n)
 * 
 * Parameters:
 * e - value to find the multiplicative inverse of
 * phi(n) - (p-1) * (q-1) where n = p*q
 * print - governs whether to produce output
 * 
 * returns - a tuple of (gcd(e,phin), d)
 */
tuple<unsigned int, unsigned int> ExtEuclidian(unsigned int e, unsigned int phin, bool print = false) {
    if (print) {
        cout << "e = " << e << endl;
        cout << "i\t|qi\t|r\t|ri+1\t|ri+2\t|si\t|ti" << endl;
    }

    int q[constant::MAX_EUCL_ITER];

    function<int(int, int*) > si = [&si] (int j, int q[constant::MAX_EUCL_ITER]) -> int { // define some lambda functions for si and ti, slows performance but ¯\_(ツ)_/¯
        if (j == 0) return 1;
        else if (j == 1) return 0;
        else {
            return si(j - 2, q) - q[j - 2] * si(j - 1, q);
        }
    };

    function<int(int, int*) > ti = [&ti] (int j, int q[constant::MAX_EUCL_ITER]) -> int {
        if (j == 0) return 0;
        else if (j == 1) return 1;
        else {
            return ti(j - 2, q) - q[j - 2] * ti(j - 1, q);
        }
    };
    ///ACTUAL ALGO:
    int ri = phin;
    int ri1 = e;
    int ri2 = 0;
    int i = 0;
    si(1, q);
    while (ri1 != 0) {

        q[i] = floor(ri / ri1);
        ri2 = ri - ri1 * q[i];

        if (print) cout << i << "\t|" << q[i] << "\t|" << ri << "\t|" << ri1 << "\t|" << ri2 << "\t|" << si(i, q) << "\t|" << ti(i, q) << endl;

        ri = ri1;
        ri1 = ri2;

        i++;
    }

    if (print) cout << i << "\t|" << q[i] << "\t|" << ri << "\t|" << ri1 << "\t|" << ri2 << "\t|" << si(i, q) << "\t|" << ti(i, q) << endl;

    int d = modulo(ti(i, q), phin);

    return make_tuple(ri, d);
}

tuple<unsigned int, unsigned int> findKeyPair(unsigned int phin, bool print = false) //error indicated by private key=-1
{
    int d = -1;
    int e = 3;
    while (d == -1 && e < phin) {
        tuple<unsigned int, unsigned int> euclidRet = ExtEuclidian(e, phin, print);

        if (get<0>(euclidRet) != 1) {
            e++;
            continue;
        }
        d = get<1>(euclidRet);
    }
    return make_tuple(e, d);
}

void printCryptSys(map<char, unsigned int> cryptMap) {
    map<char, unsigned int>::iterator it = cryptMap.begin(); // iterate through map
    while (it != cryptMap.end()) {
        cout << it->first << " = " << it->second << ", ";
        it++;
    }
    it = cryptMap.begin();
    while (it != cryptMap.end()) {
        cout << endl << it->first << " = " << bitset<32>(it->second);
        it++;
    }
    cout << endl;
}

map <char, unsigned int > createRSASys() {
    map<char, unsigned int> cryptMap;

    unsigned int p = getPrime();
    unsigned int q = 0;

    while ((q = getPrime()) == p) {
    } //make sure p and q are different

    unsigned int n = p*q;
    unsigned int phin = (p - 1) * (q - 1);

    tuple<unsigned int, unsigned int> keyPair = findKeyPair(phin);

    unsigned int e = get<0>(keyPair);
    unsigned int d = get<1>(keyPair);

    cryptMap['p'] = p;
    cryptMap['q'] = q;
    cryptMap['n'] = n;
    cryptMap['e'] = e;
    cryptMap['d'] = d;

    return cryptMap;
}

bitset<32> hashFunc(bitset<14 * 8> r) {
    bitset < r.size() > mask(0xFF);
    bitset<14 * 8> result = 0x00;
    for (size_t N = 0; N < ceil(r.size() / 8); N++) {
        result = result ^ ((r >> N * 8) & mask);
    }
    return result.to_ulong() & mask.to_ulong();
}

bitset<32> trentSig(string name, unsigned int e, unsigned int n, unsigned int d, bool print = false) {
    const int namebits = 8 * 6;
    const int nbits = 8 * 4;
    const int ebits = 8 * 4;

    string bitString = "";
    int j = name.length() - 1;
    for (int i = 0; i < namebits / 8; i++) {
        if (j >= 0) {
            bitString = bitset<8>((char) name[j]).to_string() + bitString; // fill in name from back to front
        } else {
            bitString = bitset<8>(' ').to_string() + bitString; // pad with spaces
        }
        j--;
    }

    bitset<namebits> namebin = bitset<namebits>(bitString);
    bitset<nbits> nbin = bitset<nbits>(n);
    bitset<ebits> ebin = bitset<ebits>(e);

    //pretty bad workaround because I couldn't find a good way to resize or concat bitsets
    bitset < namebits + nbits + ebits> r = bitset < namebits + nbits + ebits > (namebin.to_string() + nbin.to_string() + ebin.to_string());

    bitset<32> hr = hashFunc(r); //hash

    unsigned int s = fastExp(hr.to_ulong(), d, n); //signature

    if (print) {
        cout << "r = " << r << endl << "h(r) = " << hr << endl << "s = " << bitset<32>(s) << endl;

        cout << endl << "line:209" << endl;
        cout << "h(r) = " << hr.to_ulong() << ", s = " << s << endl;
    }

    return bitset<32>(s);
}

unsigned int randomNumFromN(unsigned int n, bool print = false) {
    bitset<32> nbits(n);

    size_t k = nbits.size();
    while (nbits[k - 1] != 1) { //skip logic to skip padding of the binary number
        k--;
    }
    k--;

    unsigned int u = 0;
    for (size_t i = k - 1; i > 0; i--) {
        if (i == k - 1) u += 0b1 << (i - 1);
        else if (i == 1) u += 0b1;
        else {
            unsigned char b = randomBit(0); //can call with zero as the parameter in this func was only used for printing purposes
            u += b << i;
        }

    }

    if (print) {
        cout << "k = " << k << ", u = " << u << endl;
        cout << endl << "line:233" << endl;
        cout << "u = " << bitset<32>(u) << endl;
    }
    return u;
}

unsigned int encryptDecrypt(unsigned int message, unsigned int e_or_d, unsigned int n, bool print = false) {
    return fastExp(message, e_or_d, n, print);
}

int main(int argc, char** argv) {
    for (int run = 1; run <= 20; run++) {
        cout<<endl<<"|delim|"<<endl; //for easy file writing
        
        cout << "line:124" << endl;
        unsigned int potPrime = randomNumber(true);
        cout << "Number|" << potPrime << "|" << bitset<32>(potPrime) << endl;

        while (isPrime(potPrime)) {
            potPrime = randomNumber(); // ensure that we have a non-prime number for printing purposes
        }

        cout << endl << "line:139" << endl; // line 139 logic
        isPrime(potPrime, true);

        cout << endl << "line:145" << endl; //line 145 logic
        potPrime = getPrime();
        isPrime(potPrime, true); // print result for a prime number 
        cout << potPrime << " is perhaps a prime" << endl;

        //end line 145 logic


        //line 162 logic
        cout << endl << "line:162" << endl;
        unsigned int p = potPrime; // p assignment
        while (potPrime == p) { // make sure p and q are different
            potPrime = getPrime();
        }
        unsigned int q = potPrime; // q assignment
        unsigned int n = p*q;
        unsigned int phin = (p - 1) * (q - 1);

        tuple<unsigned int, unsigned int> keyPair = findKeyPair(phin, true);
        while (get<1>(keyPair) == -1) {//make sure we find a public key that is relatively prime
            p = getPrime();
            while ((q = getPrime()) == p) {
            } //make sure p and q are different

            n = p*q;
            phin = (p - 1)*(q - 1);

            keyPair = findKeyPair(phin);
        }
        //end line 162 logic

        //line 173
        cout << endl << "line:173" << endl;
        cout << "d = " << get<1>(keyPair) << endl;
        //end line 173

        //line 177 logic
        unsigned int e = get<0>(keyPair);
        unsigned int d = get<1>(keyPair);

        map<char, unsigned int> cryptMap = createRSASys();

        cout << endl << "line:177" << endl;
        printCryptSys(cryptMap);
        //end line 177 logic

        //line 185 logic
        map<char, unsigned int> aliceRSA = cryptMap;
        map<char, unsigned int> trentRSA = createRSASys();

        cout << endl << "line:185" << endl;
        printCryptSys(trentRSA);
        //end 185 logic

        //line 207&209 logic
        cout << endl << "line:207" << endl;
        bitset<32> sig = trentSig("Alice", aliceRSA['e'], aliceRSA['n'], trentRSA['d'], true);
        //end line 207 logic

        //line 231&233 logic
        cout << endl << "line:231" << endl;
        unsigned int u = randomNumFromN(aliceRSA['n'], true);
        //end line 231&233 logic

        //line 239 logic
        cout << endl << "line:239" << endl;
        bitset<32> hu(hashFunc(bitset<14 * 8>(u)));

        unsigned int v = encryptDecrypt(hu.to_ulong(), aliceRSA['d'], aliceRSA['n']); //private key decryption by Alice
        unsigned int Ev = encryptDecrypt(v, aliceRSA['e'], aliceRSA['n']); //public key encryption done by Bob
        cout << "u = " << u << ", h(u) = " << hu.to_ulong() << ", v = " << v << ", Ev = " << Ev << endl;
        //end line 239 logic

        //line 242 logic
        cout << endl << "line:242" << endl;
        encryptDecrypt(v, aliceRSA['e'], aliceRSA['n'], true); //repeated for printing
        //end line 242 logic

    }
    //trentSig("Alice", 136, aliceRSA['n'], trentRSA['d'], true);
    return 0;
}

