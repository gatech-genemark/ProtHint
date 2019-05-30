#ifndef KERNEL_H
#define KERNEL_H

/// Abstract Kernel class
class Kernel {
public:
    virtual double getWeight(int offset) = 0;
    virtual void setWidth(int width);
    virtual ~Kernel() {}
protected:
    int width;
};

/// Linear Kernel
class LinearKernel : public Kernel {
public:
    double getWeight(int offset);
};

/// Triangular Kernel
class TriangularKernel : public Kernel {
public:
    double getWeight(int offset);
};

/// Parabolic Kernel
class ParabolicKernel : public Kernel {
public:
    double getWeight(int offset);
};

/// Triweight Kernel
class TriweightKernel : public Kernel {
public:
    double getWeight(int offset);
};

#endif /* KERNEL_H */
