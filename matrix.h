#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <type_traits>

size_t Max (const size_t& first, const size_t& second) { return first > second ? first : second; }

class BigInteger{

    public:

        BigInteger () : _size(1), _sign(1), _digit(1, 0) {}
        BigInteger (int);
        BigInteger (const BigInteger& old) : _size(old._size), _sign(old._sign), _digit(old._digit) {}
        BigInteger (const std::string& s);

        BigInteger& operator= (const BigInteger&);

        void incOrDec (const BigInteger& other, bool mode);
        BigInteger& operator+= (const BigInteger& other) { incOrDec(other, false); return *this; }
        BigInteger& operator-= (const BigInteger& other) { incOrDec(other, true); return *this; }

        BigInteger operator- () const;

        BigInteger& operator*= (const BigInteger& other);
        BigInteger& operator*= (int other);

        BigInteger& operator/= (BigInteger other);
        BigInteger& operator%= (const BigInteger& other);

        friend std::ostream& operator<< (std::ostream&, const BigInteger&);
        friend std::istream& operator>> (std::istream&, BigInteger&);

        std::string toString() const;

        explicit operator bool () const { return ! _isZero(); }

        bool isAbsoluteLessThen (const BigInteger&) const;

        bool isAbsoluteLessOrEqualThen (const BigInteger& other) const { return isAbsoluteLessThen(other) || _digit == other._digit; }

        bool isSignEqualTo (const BigInteger& other) const { return _sign == other._sign; }

        bool isSignLessThen (const BigInteger& other) const { return _sign < other._sign; }

        bool isNegative () const { return _sign == -1; }

        bool isDigitEqualTo (const BigInteger& other) const { return _digit == other._digit; }

    private:

        size_t _size;
        int _sign;
        std::vector<int> _digit;

        static const long long _BASE = 1000000000;
        static const size_t _ELEMENT_LENGTH = 9;

        void _setN0 (const size_t& length);

        void _sumOfAbsolute (const BigInteger& other, const int& resultSign);

        void _subtractionOfAbsolute (const BigInteger& other, const int& resultSign, const int& reverseSign);

        bool _isZero () const { return _size == 1 && _digit[0] == 0; }

        void _removeZeros ();

        void _makeCorrect();

        int binSearch(const BigInteger& other);

        std::string _wholeDigit(const int& index) const;

        BigInteger _div2() const;
};

BigInteger::BigInteger (int integer){
    _sign = integer >= 0 ? 1 : -1;
    _digit.clear();
    if (integer < 0)
        integer = -integer;
    do{
        _digit.emplace_back(integer % _BASE);
        integer /= _BASE;
    } while (integer);
    _size = _digit.size();
}

BigInteger& BigInteger::operator= (const BigInteger& old){
    if (this == &old) return *this;
    _size = old._size;
    _sign = old._sign;
    _digit = old._digit;
    return *this;
}

void BigInteger::incOrDec(const BigInteger& other, bool mode){
    if ((_sign == other._sign) ^ mode)
        _sumOfAbsolute(other, _sign);
    else{
        int resultSign = _sign;
        int reverseSign = 1;
        if (isAbsoluteLessThen(other)){
            resultSign = (mode ? -1 : 1) *  other._sign;
            reverseSign = -1;
        }
        _subtractionOfAbsolute(other, resultSign, reverseSign);
    }
}

BigInteger BigInteger::operator- () const{
    BigInteger negativeCopy = *this;
    if (!_isZero()) negativeCopy._sign = -negativeCopy._sign;
    return negativeCopy;
}

BigInteger& BigInteger::operator*= (const BigInteger& other){
    BigInteger myCopy = *this;
    _setN0(myCopy._size + other._size);

    for (size_t i = 0; i < myCopy._size; ++i){
        long long piggyBank = 0;
        long long multiplier = myCopy._digit[i];

        for (size_t j = 0; j < other._size || piggyBank; ++j){
            piggyBank += _digit[i + j];
            if (j < other._size)
                piggyBank += multiplier * other._digit[j];

            _digit[i + j] = piggyBank % _BASE;
            piggyBank /= _BASE;
        }
    }

    _sign = myCopy._sign * other._sign;
    _removeZeros();
    return *this;
}

BigInteger& BigInteger::operator*= (int other){
    int last = 0;
    if (other < 0)
        _sign *= -1;
    if (other == 0)
        _sign = 0;
    other = abs(other);
    for (size_t i = 0; i < _size; ++i){
        long long temp = (long long)_digit[i] * other + (long long)last;
        _digit[i] = temp % _BASE;
        last = temp / _BASE;
    }
    if (last != 0) {
        _digit.push_back(last);
        _size++;
    }
    return *this;
}

std::string BigInteger::toString() const{
    std::string result = "";
    if (_sign == -1)
        result += "-";
    result += std::to_string(_digit[_size - 1]);

    if (_size > 2){
        for (size_t index = _size - 2; index > 0; --index){
            result += _wholeDigit(index);
        }
    }

    if (_size > 1)
        result += _wholeDigit(0);

    return result;
}

BigInteger::BigInteger (const std::string& line){
    _digit.clear();

    size_t intBegin = 0;
    if (line[intBegin] == '-'){
        _sign = -1;
        intBegin++;
    }
    else
        _sign = 1;

    size_t currentLength = 0;
    unsigned int currentDigit = 0;
    unsigned int powerOf10 = 1;
    for (size_t lastAdded = line.length() - 1; lastAdded > intBegin; lastAdded--){
        currentDigit += powerOf10 * (line[lastAdded] - '0');
        currentLength++;
        powerOf10 *= 10;
        if (currentLength == _ELEMENT_LENGTH){
            _digit.emplace_back(currentDigit);
            currentDigit = 0;
            currentLength = 0;
            powerOf10 = 1;
        }
    }

    currentDigit += powerOf10 * (line[intBegin] - '0');
    _digit.emplace_back(currentDigit);
    _size = _digit.size();
}


bool BigInteger::isAbsoluteLessThen (const BigInteger& other) const{
    if (_size != other._size)
        return _size < other._size;

    for (size_t i = _size - 1; i > 0; --i){
        if (_digit[i] != other._digit[i])
            return _digit[i] < other._digit[i];
    }

    return _digit[0] < other._digit[0];
}

void BigInteger::_setN0 (const size_t& length){
    _size = length;
    _sign = 1;
    _digit.assign(length, 0);
}

void BigInteger::_sumOfAbsolute (const BigInteger& other, const int& resultSign){
    long long piggyBank = 0;
    _sign = resultSign;
    size_t resultSize = Max(_size, other._size) + 1;
    _digit.resize(resultSize);

    for (size_t i = 0; i < resultSize; ++i){
        if (i <_size)
            piggyBank += _digit[i];
        if (i < other._size){
            piggyBank += other._digit[i];
        }

        if (piggyBank >= _BASE){
            _digit[i] = piggyBank - _BASE;
            piggyBank = 1;
        }
        else{
            _digit[i] = piggyBank;
            piggyBank = 0;
        }
    }

    _size = resultSize;
    _removeZeros();
}

void BigInteger::_subtractionOfAbsolute (const BigInteger& other, const int& resultSign, const int& reverseSign){
    long long piggyBank = 0;
    _sign = resultSign;
    size_t resultSize = Max(_size, other._size);
    _digit.resize(resultSize);

    for (size_t i = 0; i < resultSize; ++i){
        if (i < _size)
            piggyBank += reverseSign * _digit[i];
        if (i < other._size){
            piggyBank -= reverseSign * other._digit[i];
        }

        if (piggyBank < 0){
            _digit[i] = piggyBank + _BASE;
            piggyBank = -1;
        }
        else{
            _digit[i] = piggyBank;
            piggyBank = 0;
        }
    }

    _size = resultSize;
    _removeZeros();
}

void BigInteger::_removeZeros (){
    while (_size > 1 && _digit.back() == 0){
        _digit.pop_back();
        _size--;
    }

    if (_size == 1 && _digit[0] == 0)
        _sign = 1;
}

std::string BigInteger::_wholeDigit(const int& index) const{
    std::string result = std::to_string(_digit[index]);
    while (result.length() < _ELEMENT_LENGTH)
        result = "0" + result;
    return result;
}

BigInteger BigInteger::_div2() const{
    int divisor = 2;
    BigInteger result = *this;
    long long piggyBank = 0;

    for (int i = result._size - 1; i >= 0; --i){
        piggyBank += result._digit[i];
        result._digit[i] = piggyBank / divisor;
        piggyBank = (piggyBank % divisor) * _BASE;
    }

    result._removeZeros();
    return result;
}

bool operator== (const BigInteger& first, const BigInteger& second){
    return first.isSignEqualTo(second) && first.isDigitEqualTo(second);
}

bool operator!= (const BigInteger& first, const BigInteger& second){
    return !(first == second);
}

bool operator< (const BigInteger& first, const BigInteger& second){
    if (! first.isSignEqualTo(second))
        return first.isSignLessThen(second);
    else
        return first.isNegative() ? second.isAbsoluteLessThen(first) : first.isAbsoluteLessThen(second);
}

bool operator> (const BigInteger& first, const BigInteger& second){
    return second < first;
}

bool operator<= (const BigInteger& first, const BigInteger& second){
    return !(first > second);
}

bool operator>= (const BigInteger& first, const BigInteger& second){
    return !(first < second);
}

BigInteger operator+ (const BigInteger& first, const BigInteger& second){
    BigInteger result = first;
    result += second;
    return result;
}

BigInteger& operator++ (BigInteger& current){
    current += 1;
    return current;
}
BigInteger operator++ (BigInteger& current, int){
    BigInteger old = current;
    current += 1;
    return old;
}

BigInteger operator- (const BigInteger& first, const BigInteger& second){
    BigInteger result = first;
    result -= second;
    return result;
}

BigInteger& operator-- (BigInteger& current){
    current -= 1;
    return current;
}
BigInteger operator-- (BigInteger& current, int){
    BigInteger old = current;
    current -= 1;
    return old;
}

BigInteger operator* (const BigInteger& first, const BigInteger& second){
    BigInteger result = first;
    result *= second;
    return result;
}

BigInteger operator* (const BigInteger& first, const int& second){
    BigInteger result = first;
    result *= second;
    return result;
}

void BigInteger::_makeCorrect() {
    while (_digit.size() != 1 && _digit.back() == 0)
        _digit.pop_back();
    if (_isZero())
        _sign = 1;
}

int BigInteger:: binSearch(const BigInteger& other) {
    int left = 0, right = _BASE;
    while (right > left + 1) {
        int mid = (right + left) / 2;
        if (other * mid <= *this) {
            left = mid;
        }
        else {
            right = mid;
        }
    }
    return left;
}

BigInteger& BigInteger:: operator/=(BigInteger other) {

    //std :: cout << *this << std :: endl;
    //std :: cout << other << std :: endl;
    int ansSign = _sign * other._sign;
    _sign = 1;
    other._sign = 1;
    if (*this < other) {
        *this = 0;
        return *this;
    }
    BigInteger current;
    for (size_t i = 0; i != other._digit.size() - 1; i++) {
        current = current * _BASE + _digit[_digit.size() - i - 1];
    }
    BigInteger result;
    for (size_t i = other._digit.size() - 1; i != _digit.size(); i++) {
        current *= _BASE;
        current += _digit[_digit.size() - i - 1];
        int digit = current.binSearch(other);
        current -= other * digit;
        result *= _BASE;
        result += digit;
    }
    *this = result;
    _sign = ansSign;
    _makeCorrect();
    return *this;
}

BigInteger operator/ (const BigInteger& first, const BigInteger& second){
    BigInteger result = first;
    result /= second;
    //std :: cout << first << " " << second << " " << result << std :: endl;
    return result;
}

BigInteger& BigInteger::operator%= (const BigInteger& other) {
    *this = *this - (*this / other) * other;
    return *this;
}
BigInteger operator% (const BigInteger& first, const BigInteger& second){
    BigInteger result = first;
    result %= second;
    return result;
}

std::ostream& operator<< (std::ostream& out, const BigInteger& toPrint){
    out << toPrint.toString();
    return out;
}

std::istream& operator>> (std::istream& in, BigInteger& toRead){
    std::string line;
    in >> line;
    toRead = line;
    return in;
}



void MakePositive(BigInteger& bigint){
    if (bigint < 0)
        bigint *= -1;
}

BigInteger GCD (const BigInteger& first, const BigInteger& second) {
    BigInteger _pair[] = {first, second};
    MakePositive(_pair[0]); MakePositive(_pair[1]);
    std::pair<int, int> index = {0, 1};
    while (_pair[index.first]){
        _pair[index.second] %= _pair[index.first];
        std::swap(index.first, index.second);
    }

    if (first <= 0 && second <= 0)
        _pair[index.second] *= -1;
    return _pair[index.second];
}

BigInteger BinPow(const BigInteger& bigint, unsigned int grad){
    BigInteger result = 1;
    BigInteger bigintInCurrentGrad = bigint;
    for (unsigned long long currentGrad = 1; currentGrad <= grad; currentGrad <<= 1){
        if (currentGrad & grad){
            result *= bigintInCurrentGrad;
        }
        bigintInCurrentGrad = bigintInCurrentGrad * bigintInCurrentGrad;
    }
    return result;
}

double ToDouble(const std::string& stringDouble = "0"){
    double beforePoint = 0;
    double afterPoint = 0;
    size_t index = 0;
    double sign = 1;
    if (stringDouble[0] == '-'){
        sign = -1;
        index = 1;
    }
    while (index < stringDouble.size() && stringDouble[index] != '.'){
        double digit = stringDouble[index] - '0';
        beforePoint = beforePoint * 10 + digit;
        ++index;
    }

    size_t dotIndex = index;
    index = stringDouble.size() - 1;
    while (index > dotIndex){
        double digit = stringDouble[index] - '0';
        afterPoint = (afterPoint + digit) * 0.1;
        --index;
    }

    return sign * (beforePoint + afterPoint);
}

class Rational{

    public:

        Rational (int integer = 0) : _numerator(integer), _denominator(1) {}

        Rational (const BigInteger& bigint) : _numerator(bigint), _denominator(1) {}

        Rational (const Rational& old) : _numerator(old._numerator), _denominator(old._denominator) {}

        Rational& operator= (const Rational& old) {
            _numerator = old._numerator;
            _denominator = old._denominator;
            return *this;
        }

        Rational& operator+= (const Rational& other) {
            _numerator *= other._denominator;
            _numerator += other._numerator * _denominator;
            _denominator *= other._denominator;
            Normalize();
            return *this;
        }

        Rational& operator-= (const Rational& other) {
            _numerator *= other._denominator;
            _numerator -= other._numerator * _denominator;
            _denominator *= other._denominator;
            Normalize();
            //std :: cout << _numerator._size << " " << _denominator._size << std :: endl;
            return *this;
        }

        Rational& operator*= (const Rational& other) {
            _numerator *= other._numerator;
            _denominator *= other._denominator;
            Normalize();
            return *this;
        }

        Rational operator- () {
            Rational result = *this;
            result *= -1;
            return result;
        }

        Rational& operator/= (const Rational& other) {
            _numerator *= other._denominator;
            _denominator *= other._numerator;
            Normalize();
            return *this;
        }

        bool isLess (const Rational& other) const {
            return _numerator * other._denominator < other._numerator * _denominator;
        }

        bool isEqual (const Rational& other) const {
            return _numerator == other._numerator && _denominator == other._denominator;
        }

        void Normalize() {
            MakeCoprime();
            if (_denominator < 0){
                _numerator *= -1;
                _denominator *= -1;
            }
        }

        void MakeCoprime() {
            BigInteger gcd = GCD(_numerator, _denominator);
            _numerator /= gcd;
            _denominator /= gcd;
        }

        std::string toString() const{
            std::string result = _numerator.toString();
            if (_denominator != 1)
                result += "/" + _denominator.toString();
            return result;
        }

        std::string asDecimal(size_t precision = 0) const {
            const BigInteger base = 10;
            BigInteger intResult = (_numerator * BinPow(base, precision)) / _denominator;
            std::string result = intResult.toString();

            std::string sign = "";
            if (result[0] == '-'){
                sign = "-";
                result = result.substr(1, result.size() - 1);
            }
            while (result.size() < precision + 1)
                result = "0" + result;

            size_t decimalPart = result.size() - precision;
            if (precision)
                return sign + result.substr(0, decimalPart) + "." + result.substr(decimalPart, precision);
            else
                return sign + result;
        }

        explicit operator double() const{
            const size_t precision = 15;
            return ToDouble(asDecimal(precision));
        }

    //private:
        BigInteger _numerator;
        BigInteger _denominator;
};

Rational operator+ (const Rational& first, const Rational& second){
    Rational result = first;
    result += second;
    return result;
}

Rational operator- (const Rational& first, const Rational& second){
    Rational result = first;
    result -= second;
    return result;
}

Rational operator* (const Rational& first, const Rational& second){
    Rational result = first;
    result *= second;
    return result;
}

Rational operator/ (const Rational& first, const Rational& second){
    Rational result = first;
    result /= second;
    return result;
}

std::istream& operator>> (std::istream& in, Rational& toRead){
    std::string line;
    in >> line;
    std :: string numeratorStr;
    BigInteger denumerator = 1;
    for (unsigned i = 0; i < line.size(); i++) {
        if (line[i] != '.')
            numeratorStr.push_back(line[i]);
        else
            denumerator = BinPow(10, line.size() - i - 1);
    }
    BigInteger numerator = numeratorStr;
    toRead = (Rational)(numerator) / (Rational)(denumerator);
    return in;
}

bool operator== (const Rational& first, const Rational& second){
    return first.isEqual(second);
}

bool operator!= (const Rational& first, const Rational& second){
    return !(first == second);
}

bool operator< (const Rational& first, const Rational& second){
    return first.isLess(second);
}

bool operator> (const Rational& first, const Rational& second){
    return second < first;
}

bool operator<= (const Rational& first, const Rational& second){
    return !(first > second);
}

bool operator>= (const Rational& first, const Rational& second){
    return !(first < second);
}


template <unsigned N>
struct mySqrt {
    static const unsigned value = sqrt(N);
};

template <unsigned N, unsigned div>
struct IsPrimeRec {
    static const bool value = (N % div != 0) && IsPrimeRec <N, div - 1>::value;
};

template <unsigned N>
struct IsPrimeRec <N, 1> {
    static const bool value = true;
};

template <unsigned N>
struct IsPrime {
    static const bool value = IsPrimeRec <N, mySqrt <N> :: value>::value;
};

template <>
struct IsPrime <1> {
    static const bool value = false;
};

struct CompileError {
    virtual void f();
};

template <bool arg>
struct Check {};

template <>
struct Check<false> {
    CompileError cmpEr;
};

template <unsigned N>
class Finite {
public:
    Finite(){}

    Finite(int x);

    Finite(const Finite& other);

    Finite& operator =(const Finite& other);

    Finite& operator +=(const Finite& other);

    Finite& operator ++();

    Finite& operator -=(const Finite& other);

    Finite& operator *=(const Finite& other);

    Finite pow(unsigned k) const;

    Finite& operator /=(const Finite& other);

    Finite inverse() const;

    bool operator !=(const Finite& other) const;

    void print() const{
        std :: cout << _x;
    }
private:
    unsigned _x;
};

template <unsigned N>
Finite<N>:: Finite(int x) {
    x %= (signed)N;
    if (x < 0)
        x += N;
    _x = x;
}

template <unsigned N>
Finite<N>:: Finite(const Finite& other) {
    //std :: cout << other._x << std :: endl;
    _x = other._x;
}

template <unsigned N>
Finite<N>& Finite<N> :: operator =(const Finite& other) {
    _x = other._x;
    return *this;
}

template <unsigned N>
Finite<N>& Finite<N> :: operator +=(const Finite& other) {
    _x += other._x;
    if (_x >= N)
        _x -= N;
    return *this;
}

template <unsigned N>
Finite<N>& Finite<N> :: operator ++() {
    *this += 1;
    return *this;
}

template <unsigned N>
Finite<N> operator +(const Finite<N>& first, const Finite<N>& second) {
    Finite<N> result = first;
    result += second;
    return result;
}

template <unsigned N>
Finite<N>& Finite<N> :: operator -=(const Finite& other) {
    _x += N - other._x;
    if (_x >= N)
        _x -= N;
    return *this;
}

template <unsigned N>
Finite<N> operator -(const Finite<N>& first, const Finite<N>& second) {
    Finite<N> result = first;
    result -= second;
    return result;
}

template <unsigned N>
Finite<N>& Finite<N> :: operator *=(const Finite& other) {
    _x = (_x * other._x) % N;
    //std :: cout << _x << std :: endl;
    return *this;
}

template <unsigned N>
Finite<N> operator *(const Finite<N>& first, const Finite<N>& second) {
    Finite<N> result = first;
    result *= second;
    return result;
}

template <unsigned N>
Finite<N> Finite<N>:: pow(unsigned k) const {
    Finite result = 1;
    Finite mult = *this;
    while (k > 0) {
        if (k & 1) {
            result *= mult;
        }
        mult *= mult;
        k >>= 1;
    }
    return result;
}

template <unsigned N>
Finite<N> Finite<N>:: inverse() const {
    Check<IsPrime<N> :: value> ch __attribute__((unused));
    return pow(N - 2);
}

template <unsigned N>
Finite<N>& Finite<N> :: operator /=(const Finite& other) {
    _x = (_x * other.inverse()._x) % N;
    //std :: cout << _x;
    return *this;
}

template <unsigned N>
Finite<N> operator /(const Finite<N>& first, const Finite<N>& second) {
    Finite<N> result = first;
    result /= second;
    return result;
}

template <unsigned N>
bool Finite<N> :: operator !=(const Finite<N>& other) const {
    return _x != other._x;
}



template <typename Field>
struct IsFieldCorrect {};

template <unsigned N>
struct IsFieldCorrect <Finite<N>> {
    Check<IsPrime<N> :: value> ch;
};

template <unsigned N, unsigned M>
struct CheckSquare {
    CompileError cmpEr;
};

template <unsigned N>
struct CheckSquare<N, N> {};

template <typename Field>
class Wrapper{
public:
    Wrapper(std::vector <Field>* _a) {
        a = _a;
    }

    Field& operator[](const unsigned &ind) {
        return (*a)[ind];
    }

private:
    std :: vector <Field>* a;
};

template <typename Field>
class ConstWrapper{
public:
    ConstWrapper(const std::vector <Field>* _a) {
        a = _a;
    }

    const Field& operator[](const unsigned &ind) {
        return (*a)[ind];
    }
private:
    const std :: vector <Field>* a;
};

template <unsigned M, unsigned N, typename Field = Rational>
class Matrix : private IsFieldCorrect<Field> {
public:
    Matrix();

    Matrix(std::vector <std :: vector <Field> > a);

    Matrix(std::initializer_list<std::initializer_list<int>> a);

    Matrix(const Matrix& other);

    bool operator != (const Matrix& other) const;

    bool operator == (const Matrix& other) const;

    Matrix& operator =(const Matrix& other);

    Matrix& operator +=(const Matrix& other);

    Matrix& operator -=(const Matrix& other);

    Matrix& operator *=(const Field& scal);

    Matrix& operator *=(const Matrix<N, N, Field>& other);

    Field det() const;

    Matrix<N, M, Field> transposed() const;

    unsigned rank() const;

    Matrix& invert();

    Matrix inverted() const;

    Field trace() const;

    std :: vector <Field> getRow(unsigned k) const;

    std :: vector <Field> getColumn(unsigned k) const;

    Wrapper<Field> operator[](const unsigned& ind);

    ConstWrapper<Field> operator[](const unsigned& ind) const;

    std :: vector <std :: vector <Field> > makeCorrectSize(unsigned newN, unsigned newM) const {
        std::vector<std::vector<Field> > result = _a;

        while (result.size() < newN) {
            result.push_back({});
        }
        while (result.size() > newN) {
            result.pop_back();
        }

        for (unsigned i = 0; i < newN; ++i) {
            result[i].resize(newM);
        }

        return result;
    }

private:

    std::vector <std :: vector <Field> > _a;
    IsFieldCorrect<Field> go;

    void directGauss(Field& sign, unsigned lim = N) {
        for (unsigned i = 0; i < std :: min(M, N); i++) {
            for (unsigned j = i + 1; j < M; j++) {
                if (_a[j][i] != 0) {
                    std :: swap(_a[j], _a[i]);
                    sign *= -1;
                    break;
                }
            }
            //int c = 0;
            //std :: cout << i << std :: endl;
            if (_a[i][i] != 0) {
                for (unsigned j = i + 1; j < M; j++) {
                    Field now = _a[j][i] / _a[i][i];
                    for (unsigned k = i; k < lim; k++) {
                        _a[j][k] -= now * _a[i][k];
                        /*if (abs(_a[j][k]) < 1e-6)
                            _a[j][k] = 0;*/
                        //std :: cout << c++ << std :: endl;
                    }
                }
            }
            /*for (unsigned j = 0; j < M; j++) {
                for (unsigned k = 0; k < lim; k++) {
                    std :: cout << _a[j][k] << " ";
                }
                std :: cout << std :: endl;
            }
            std :: cout << "!!!!!!!!!" << std :: endl;*/
        }
    }

    void inverseGauss() {
        for (unsigned i = M - 1; i + 1 != 0; i--) {
            //std :: cout << i << std :: endl;
            Field now = _a[i][i];
            for (unsigned j = i; j < 2 * N; j++) {
                _a[i][j] /= now;
            }
            for (unsigned j = i - 1; j + 1 != 0; j--) {
                Field now = _a[j][i];
                for (unsigned k = i; k < 2 * N; k++) {
                    //std :: cout << i << " " << j << " " << k << std :: endl;
                    _a[j][k] -= now * _a[i][k];
                    /*if (abs(_a[j][k]) < 1e-6)
                            _a[j][k] = 0;*/
                }
            }
        }
    }
};

template <unsigned M, unsigned N, typename Field>
std :: vector <Field> Matrix <M, N, Field> :: getRow(unsigned k) const {
    return _a[k];
}

template <unsigned M, unsigned N, typename Field>
std :: vector <Field> Matrix <M, N, Field> :: getColumn(unsigned k) const {
    std :: vector <Field> res;
    for (unsigned i = 0; i < M; i++)
        res.push_back(_a[i][k]);
    return res;
}

template <unsigned M, unsigned N, typename Field>
Matrix <M, N, Field>:: Matrix() {
    _a.resize(M, std :: vector <Field>(N, 0));
    for (unsigned i = 0; i < std :: min(M, N); i++)
        _a[i][i] = 1;
}

template <unsigned M, unsigned N, typename Field>
Matrix <M, N, Field>:: Matrix(std::vector <std :: vector <Field> > a) {
    _a.resize(M, std :: vector <Field>(N));
    for (unsigned i = 0; i < M; i++) {
        for (unsigned j = 0; j < N; j++) {
            _a[i][j] = a[i][j];
        }
    }
}

template <unsigned M, unsigned N, typename Field>
Matrix <M, N, Field>:: Matrix(std::initializer_list<std::initializer_list<int>> a) {
    _a.resize(M, std :: vector <Field>(N));
    unsigned i = 0, j = 0;
    for (std::initializer_list<std::initializer_list<int>>::iterator it1 = a.begin(); it1 != a.end(); it1++) {
        j = 0;
        for (std::initializer_list<int>::iterator it2 = (*it1).begin(); it2 != (*it1).end(); it2++) {
            //std :: cout << i << " " << j << " " << *it2 << " " << (*it1).size() << std :: endl;
            _a[i][j] = *it2;
            j++;
        }
        i++;
    }
    //std :: cout << "????" << std :: endl;
}

template <unsigned M, unsigned N, typename Field>
Matrix <M, N, Field>:: Matrix(const Matrix<M, N, Field>& other) {
    _a.resize(M, std :: vector <Field>(N));
    for (unsigned i = 0; i < M; i++) {
        for (unsigned j = 0; j < N; j++) {
            _a[i][j] = other._a[i][j];
        }
    }
}

template <unsigned M, unsigned N, typename Field>
Wrapper<Field> Matrix <M, N, Field>::operator[](const unsigned& ind) {
    return &_a[ind];
}

template <unsigned M, unsigned N, typename Field>
ConstWrapper<Field> Matrix <M, N, Field>::operator[](const unsigned& ind) const {
    return &_a[ind];
}

template <unsigned M, unsigned N, typename Field>
bool Matrix <M, N, Field> :: operator ==(const Matrix<M, N, Field>& other) const {
    for (unsigned i = 0; i < M; i++) {
        for (unsigned j = 0; j < N; j++) {
            if (_a[i][j] != other[i][j])
                return false;
        }
    }
    return true;
}

template <unsigned M, unsigned N, typename Field>
bool Matrix <M, N, Field> :: operator !=(const Matrix<M, N, Field>& other) const {
    return !(*this == other);
}

template <unsigned M, unsigned N, typename Field>
Matrix <M, N, Field>& Matrix <M, N, Field> :: operator =(const Matrix<M, N, Field>& other) {
    for (unsigned i = 0; i < M; i++) {
        for (unsigned j = 0; j < N; j++) {
            _a[i][j] = other._a[i][j];
        }
    }
    return *this;
}

template <unsigned M, unsigned N, typename Field>
Matrix <M, N, Field>& Matrix <M, N, Field> :: operator +=(const Matrix<M, N, Field>& other) {
    for (unsigned i = 0; i < M; i++) {
        for (unsigned j = 0; j < N; j++) {
            _a[i][j] += other._a[i][j];
        }
    }
    return *this;
}

template <unsigned M, unsigned N, typename Field>
Matrix <M, N, Field> operator +(const Matrix<M, N, Field>& first, const Matrix<M, N, Field>& second) {
    Matrix<M, N, Field> result = first;
    result += second;
    return result;
}

template <unsigned M, unsigned N, typename Field>
Matrix <M, N, Field>& Matrix <M, N, Field> :: operator -=(const Matrix<M, N, Field>& other) {
    for (unsigned i = 0; i < M; i++) {
        for (unsigned j = 0; j < N; j++) {
            _a[i][j] -= other._a[i][j];
        }
    }
    return *this;
}

template <unsigned M, unsigned N, typename Field>
Matrix <M, N, Field> operator -(const Matrix<M, N, Field>& first, const Matrix<M, N, Field>& second) {
    Matrix<M, N, Field> result = first;
    result -= second;
    return result;
}

template <unsigned M, unsigned N, typename Field>
Matrix <M, N, Field>& Matrix <M, N, Field> :: operator *=(const Field& scal) {
    for (unsigned i = 0; i < M; i++) {
        for (unsigned j = 0; j < N; j++) {
            _a[i][j] *= scal;
        }
    }
    return *this;
}

template <unsigned M, unsigned N, typename Field>
Matrix <M, N, Field> operator *(const Field& scal, const Matrix<M, N, Field>& second) {
    Matrix<M, N, Field> result = second;
    result *= scal;
    return result;
}

template <unsigned M, unsigned N, typename Field>
Matrix <M, N, Field> operator *(const Matrix<M, N, Field>& first, const Field& scal) {
    Matrix<M, N, Field> result = first;
    result *= scal;
    return result;
}

template <unsigned N, unsigned M, typename Field>
Field Matrix <N, M, Field>:: det() const {
    CheckSquare<N, M>* cS = new CheckSquare<N, M>;
    delete cS;
    Field sign = 1;
    Matrix temp = *this;
    temp.directGauss(sign);
    Field res = 1;
    for (unsigned i = 0; i < N; i++)
        res *= temp._a[i][i];
    return res * sign;
}

template <unsigned N, unsigned M, typename Field>
unsigned Matrix <N, M, Field>:: rank() const {
    //std :: cout << "???" << std :: endl;
    Matrix temp = *this;
    Field sign = 1;
    temp.directGauss(sign);
    unsigned res = 0;
    for (unsigned i = 0; i < N; i++) {
        for (unsigned j = 0; j < M; j++) {
            if (temp._a[i][j] != 0) {
                res++;
                break;
            }
        }
    }
    return res;
}

template <unsigned N, unsigned M, typename Field>
Matrix <N, M, Field>& Matrix <N, M, Field>:: invert() {
    CheckSquare<N, M>* cS = new CheckSquare<N, M>;
    delete cS;
    for (unsigned i = 0; i < N; i++) {
        for (unsigned j = 0; j < N; j++) {
            if (i != j)
                _a[i].push_back(0);
            else
                _a[i].push_back(1);
        }
    }
    Field sign = 1;
    directGauss(sign, 2 * N);
    //std :: cout << "!!!!";
    inverseGauss();
    for (unsigned i = 0; i < N; i++) {
        for (unsigned j = 0; j < N; j++) {
            _a[i][j] = _a[i][j + N];
        }
        _a[i].erase(_a[i].begin() + N, _a[i].end());
    }
    /*for (unsigned i = 0; i < N; i++) {
        for (unsigned j = 0; j < N; j++) {
            std :: cout << (double)_a[i][j];
            std :: cout << " ";
        }
        std :: cout << std :: endl;
    }*/
    return *this;
}

template <unsigned N, unsigned M, typename Field>
Matrix <N, M, Field> Matrix <N, M, Field>:: inverted() const {
    Matrix result = *this;
    result.invert();
    return result;
}

template <unsigned N, unsigned M, typename Field>
Field Matrix <N, M, Field>:: trace() const {
    CheckSquare<N, M>* cS = new CheckSquare<N, M>;
    delete cS;
    Field res = 0;
    for (unsigned i = 0; i < N; i++)
        res += _a[i][i];
    return res;
}

template <unsigned N, typename Field>
Matrix<N, N, Field> doStrassen(const Matrix<N, N, Field>& first, const Matrix<N, N, Field>& second) {
    if (N < 16) {
        Matrix<N, N, Field> result;
        for (unsigned i = 0; i < N; i++) {
            for (unsigned j = 0; j < N; j++) {
                result[i][j] = 0;
                for (unsigned k = 0; k < N; k++) {
                    result[i][j] += first[i][k] * second[k][j];
                }
            }
        }
        return result;
    }

    Matrix<N / 2, N / 2, Field> a11, a12, a21, a22, b11, b12, b21, b22, c11, c12, c21, c22, p1, p2, p3, p4, p5, p6, p7;

    for (unsigned i = 0; i < N / 2; i++) {
        for (unsigned j = 0; j < N / 2; j++) {
            a11[i][j] = first[i][j];
            a12[i][j] = first[i][j + N / 2];
            a21[i][j] = first[i + N / 2][j];
            a22[i][j] = first[i + N / 2][j + N / 2];
            b11[i][j] = second[i][j];
            b12[i][j] = second[i][j + N / 2];
            b21[i][j] = second[i + N / 2][j];
            b22[i][j] = second[i + N / 2][j + N / 2];
        }
    }

    p1 = doStrassen(a11 + a22, b11 + b22);
    p2 = doStrassen(a21 + a22, b11);
    p3 = doStrassen(a11, b12 - b22);
    p4 = doStrassen(a22, b21 - b11);
    p5 = doStrassen(a11 + a12, b22);
    p6 = doStrassen(a21 - a11, b11 + b12);
    p7 = doStrassen(a12 - a22, b21 + b22);

    c11 = p1 + p4 - p5 + p7;
    c12 = p3 + p5;
    c21 = p2 + p4;
    c22 = p1 - p2 + p3 + p6;

    Matrix<N, N, Field> result;
    for (unsigned i = 0; i < N / 2; i++) {
        for (unsigned j = 0; j < N / 2; j++) {
            result[i][j] = c11[i][j];
            result[i][j + N / 2] = c12[i][j];
            result[i + N / 2][j] = c21[i][j];
            result[i + N / 2][j + N / 2] = c22[i][j];
        }
    }
    //result[0][0].print();
    //std :: cout << std :: endl;
    return result;
}

template <unsigned N, unsigned M, unsigned K, typename Field>
Matrix<N, K, Field> operator *(const Matrix<N, M, Field>& first, const Matrix<M, K, Field>& second) {
    const unsigned minPowerOfTwo = 1 << static_cast<int>(ceil(log2(std::max(N, std::max(M, K)))));

    Matrix<minPowerOfTwo, minPowerOfTwo, Field> left = first.makeCorrectSize(minPowerOfTwo, minPowerOfTwo);
    Matrix<minPowerOfTwo, minPowerOfTwo, Field> right = second.makeCorrectSize(minPowerOfTwo, minPowerOfTwo);

    Matrix<minPowerOfTwo, minPowerOfTwo, Field> tmp = doStrassen(left, right);
    Matrix<N, K, Field> result = tmp.makeCorrectSize(N, K);

    return result;
}

template <unsigned M, unsigned N, typename Field>
Matrix <M, N, Field>& Matrix <M, N, Field> :: operator *=(const Matrix<N, N, Field>& other) {
    *this = *this * other;
    return *this;
}

template <unsigned N, unsigned M, typename Field>
Matrix <M, N, Field>  Matrix <N, M, Field>:: transposed() const {
    Matrix <M, N, Field> result;
    for (unsigned i = 0; i < N; i++) {
        for (unsigned j = 0; j < M; j++) {
            result[j][i] = _a[i][j];
            /*if (Field == Finite)
                result[j][i].print();*/
            //std :: cout << " ";
        }
        //std :: cout << std :: endl;
    }
    return result;
}

template <unsigned N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;
