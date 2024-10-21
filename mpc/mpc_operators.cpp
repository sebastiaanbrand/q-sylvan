#include <iostream>
#include <mpc.h>  // GNU MPC library

class MPC {
private:
    mpc_t value;  // multi-precision complex number

public:
    // Constructors
    MPC() {
        mpc_init2(value, 256);  // Initialize with 256 bits of precision
    }

    MPC(const char* real, const char* imag) {
        mpc_init2(value, 256);
        mpc_set_str(value, real, imag, 10);  // Initialize from string (base 10)
    }

    MPC(double real, double imag) {
        mpc_init2(value, 256);
        mpc_set_d_d(value, real, imag, MPC_RNDNN);  // Initialize from double
    }

    MPC(const MPC& other) {
        mpc_init2(value, 256);
        mpc_set(value, other.value, MPC_RNDNN);  // Copy constructor
    }

    ~MPC() {
        mpc_clear(value);  // Free the memory
    }

    // Operator overloading: Addition (+)
    MPC operator+(const MPC& other) const {
        MPC result;
        mpc_add(result.value, this->value, other.value, MPC_RNDNN);  // Perform addition
        return result;
    }

    // Operator overloading: Subtraction (-)
    MPC operator-(const MPC& other) const {
        MPC result;
        mpc_sub(result.value, this->value, other.value, MPC_RNDNN);  // Perform subtraction
        return result;
    }

    // Operator overloading: Multiplication (*)
    MPC operator*(const MPC& other) const {
        MPC result;
        mpc_mul(result.value, this->value, other.value, MPC_RNDNN);  // Perform multiplication
        return result;
    }

    // Operator overloading: Division (/)
    MPC operator/(const MPC& other) const {
        MPC result;
        mpc_div(result.value, this->value, other.value, MPC_RNDNN);  // Perform division
        return result;
    }

    // Operator overloading: Assignment (=)
    MPC& operator=(const MPC& other) {
        if (this != &other) {
            mpc_set(value, other.value, MPC_RNDNN);  // Set the value
        }
        return *this;
    }

    // Operator overloading: Equality (==)
    bool operator==(const MPC& other) const {
        return mpc_cmp(value, other.value) == 0;  // Compare two complex numbers
    }

    // Operator overloading: Inequality (!=)
    bool operator!=(const MPC& other) const {
        return !(*this == other);
    }

    // Print complex number
    void print() const {
        mpc_out_str(stdout, 10, 0, value, MPC_RNDNN);
        std::cout << std::endl;
    }

    // Access the underlying mpc_t (for advanced usage)
    const mpc_t& getValue() const {
        return value;
    }
};

int main() {
    // Create multi-precision complex numbers
    MPC z1(2.0, 3.0);       // z1 = 2 + 3i
    MPC z2(4.0, -1.0);      // z2 = 4 - 1i

    // Addition
    MPC sum = z1 + z2;
    std::cout << "z1 + z2 = ";
    sum.print();

    // Subtraction
    MPC diff = z1 - z2;
    std::cout << "z1 - z2 = ";
    diff.print();

    // Multiplication
    MPC prod = z1 * z2;
    std::cout << "z1 * z2 = ";
    prod.print();

    // Division
    MPC quot = z1 / z2;
    std::cout << "z1 / z2 = ";
    quot.print();

    // Equality and Inequality
    bool isEqual = (z1 == z2);
    std::cout << "z1 == z2: " << std::boolalpha << isEqual << std::endl;

    bool isNotEqual = (z1 != z2);
    std::cout << "z1 != z2: " << std::boolalpha << isNotEqual << std::endl;

    return 0;
}
