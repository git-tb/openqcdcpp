#ifndef PROTECTEDOBJECT_H
#define PROTECTEDOBJECT_H

#include <cstdlib>		/// std::size_t
#include <algorithm>	/// std::copy

/**
 * @brief ProtectedObjects are those that should not be changed everywhere in the code.
 * On one hand, we could just declare constant global variables at the beginning of
 * each program. Alternatively we can keep them non-const but restrict assignment by 
 * making them private in a dedicated class. This simplifies static code analysis.
 */
template <typename T>
class ProtectedObject {
private:
    T value;

public:
    ProtectedObject() = delete;
    explicit ProtectedObject(T v) : value(v) {}

    /// forbid assignment from outside
    void operator=(T) = delete;

    /// controlled setter
    void set(T v) { value = v; }

    /// implicit conversion for reading
    operator T() const { return value; }
};

template <typename T, std::size_t N>
class ProtectedObject<T[N]> {
private:
    T value[N];

public:
    ProtectedObject() = delete;
    explicit ProtectedObject(const T (&v)[N]) {
        std::copy(v, v + N, value);
    }

    void set(const T (&v)[N]) {
        std::copy(v, v + N, value);
    }

	void set(std::size_t i, T v) {
        assert(i<N);
		value[i] = v;
    }

    operator const T*() const { return value; }
};

#endif