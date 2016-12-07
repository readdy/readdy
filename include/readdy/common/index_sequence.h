/**
 * << detailed description >>
 *
 * @file index_sequence.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.08.16
 */

#ifndef READDY_MAIN_INDEX_SEQUENCE_H
#define READDY_MAIN_INDEX_SEQUENCE_H

#include <cstddef>

namespace std {
    template<size_t... Ints>
    struct index_sequence {
        using type = index_sequence;
        using value_type = size_t;

        static constexpr std::size_t size() noexcept { return sizeof...(Ints); }
    };

    // --------------------------------------------------------------

    template<class Sequence1, class Sequence2>
    struct _merge_and_renumber;

    template<size_t... I1, size_t... I2>
    struct _merge_and_renumber<index_sequence<I1...>, index_sequence<I2...>>
            : index_sequence<I1..., (sizeof...(I1) + I2)...> {
    };

    // --------------------------------------------------------------

    template<size_t N>
    struct make_index_sequence
            : _merge_and_renumber<typename make_index_sequence<N / 2>::type,
                    typename make_index_sequence<N - N / 2>::type> {
    };

    template<>
    struct make_index_sequence<0> : index_sequence<> { };
    template<>
    struct make_index_sequence<1> : index_sequence<0> { };
}
#endif //READDY_MAIN_INDEX_SEQUENCE_H
