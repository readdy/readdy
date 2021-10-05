//
// Created by mho on 11/20/19.
//

#include <catch2/catch.hpp>
#include <graphs/IndexPersistentVector.h>
#include <set>

class A {
public:

    explicit A(int x) : x(x) {}

    void deactivate() {
        active = false;
    }

    [[nodiscard]] bool deactivated() const { return !active; }

    [[nodiscard]] auto val() const { return x; }

    bool operator==(const A &other) const {
        return other.active == active && other.x == x;
    }

    bool operator<(const A &other) const {
        return other.x < x;
    }

private:
    bool active {true};
    int x;
};

bool randomBool() {
    static auto gen = std::bind(std::uniform_int_distribution<>(0,1),std::default_random_engine());
    return gen();
}

TEST_CASE("IPV usage test", "[ipv]") {
    std::vector<A> elems;
    graphs::IndexPersistentVector<A> v;

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(1, 600000);

    for(int i = 0; i < 10000; ++i) {
        REQUIRE(elems.size() == v.size());
        if(randomBool()) {
            // add a random element
            auto n = dist(rng);
            v.push_back(A(n));
            elems.emplace_back(n);
        } else {
            if(!elems.empty()) {
                std::uniform_int_distribution<std::mt19937::result_type> draw(0, elems.size()-1);
                auto it = v.begin() + draw(rng);
                auto val = it->val();
                v.erase(it);
                auto findit = std::find_if(elems.begin(), elems.end(), [val](auto a) { return a.val() == val; });
                REQUIRE(findit != elems.end());
                elems.erase(findit);
            }
        }
        REQUIRE(v.size() == elems.size());
        for(auto a : elems) {
            auto it = std::find(v.begin(), v.end(), a);
            REQUIRE(it != v.end());
            REQUIRE(it == (v.begin() + std::distance(v.begin(), it)));
            REQUIRE(it == (v.end() - std::distance(it, v.end())));
        }
        for(auto& vertex : v) {
            REQUIRE(std::find(elems.begin(), elems.end(), vertex) != elems.end());
        }
    }
}

SCENARIO("Test ipv active iterator", "[ipv]") {
    GIVEN("A IPV with five elements") {
        graphs::IndexPersistentVector<A> v;
        v.push_back(A(5));
        v.push_back(A(1));
        v.push_back(A(7));
        v.push_back(A(8));
        v.push_back(A(3));
        THEN("The distance between active begin and active end is 5") {
            REQUIRE(std::distance(v.begin(), v.end()) == 5);
        }

        WHEN("Removing element number 2") {
            v.erase(v.begin() + 2);
            THEN("The active size is 4") {
                REQUIRE(std::distance(v.begin(), v.end()) == 4);
                REQUIRE(v.size() == 4);
            }
            AND_THEN("the active iterator skips the element with number 7") {
                WHEN("Accessing it through random access") {
                    std::vector<std::size_t> mapping{0, 1, 3, 4};
                    auto ix = 0;
                    for (auto it = v.begin(); it != v.end(); ++it, ++ix) {
                        REQUIRE(it->val() == v.at(ix).val());
                    }
                }
                WHEN("Accessing it through plus operator") {
                    REQUIRE((*v.begin()).val() == 5);
                    REQUIRE((*(v.begin() + 1)).val() == 1);
                    REQUIRE((*(v.begin() + 2)).val() == 8);
                    REQUIRE((*(v.begin() + 3)).val() == 3);
                }
                WHEN("Accessing it through minus operator") {
                    REQUIRE((v.end() - 1)->val() == 3);
                    REQUIRE((v.end() - 2)->val() == 8);
                    REQUIRE((v.end() - 3)->val() == 1);
                    REQUIRE((v.end() - 4)->val() == 5);
                }
                WHEN("Accessing it through stepwise increase") {
                    auto itBegin = v.begin();
                    REQUIRE((itBegin++)->val() == 5);
                    REQUIRE((itBegin++)->val() == 1);
                    REQUIRE((itBegin++)->val() == 8);
                    REQUIRE((itBegin++)->val() == 3);
                    REQUIRE(itBegin == v.end());
                }
                WHEN("Accessing it through stepwise decrease") {
                    auto itEnd = v.end();
                    REQUIRE((--itEnd)->val() == 3);
                    REQUIRE((--itEnd)->val() == 8);
                    REQUIRE((--itEnd)->val() == 1);
                    REQUIRE((--itEnd)->val() == 5);
                    REQUIRE(itEnd == v.begin());
                }
            }
            AND_WHEN("Removing element number 1") {
                v.erase(v.begin() + 1);
                THEN("The active size is 3") {
                    REQUIRE(std::distance(v.begin(), v.end()) == 3);
                    REQUIRE(v.size() == 3);
                }
                AND_THEN("the active iterator skips the elements with numbers 7, 1") {
                    WHEN("Accessing it through random access") {
                        std::vector<std::size_t> mapping{0, 3, 4};
                        auto ix = 0;
                        for (auto it = v.begin(); it != v.end(); ++it, ++ix) {
                            REQUIRE(it->val() == v.at(ix).val());
                        }
                    }
                    WHEN("Accessing it through plus operator") {
                        REQUIRE((*v.begin()).val() == 5);
                        REQUIRE((*(v.begin() + 1)).val() == 8);
                        REQUIRE((*(v.begin() + 2)).val() == 3);
                    }
                    WHEN("Accessing it through minus operator") {
                        REQUIRE((v.end() - 1)->val() == 3);
                        REQUIRE((v.end() - 2)->val() == 8);
                        REQUIRE((v.end() - 3)->val() == 5);
                    }
                    WHEN("Accessing it through stepwise increase") {
                        auto it = v.begin();
                        REQUIRE((it++)->val() == 5);
                        REQUIRE((it++)->val() == 8);
                        REQUIRE((it++)->val() == 3);
                        REQUIRE(it == v.end());
                    }
                    WHEN("Accessing it through stepwise decrease") {
                        auto itEnd = v.end();
                        REQUIRE((--itEnd)->val() == 3);
                        REQUIRE((--itEnd)->val() == 8);
                        REQUIRE((--itEnd)->val() == 5);
                        REQUIRE(itEnd == v.begin());
                    }
                }

                AND_WHEN("Removing element number 0") {
                    v.erase(v.begin());
                    THEN("The active size is 2") {
                        REQUIRE(std::distance(v.begin(), v.end()) == 2);
                        REQUIRE(v.size() == 2);
                    }

                    AND_THEN("the active iterator skips the elements with numbers 7, 1, 5") {
                        WHEN("Accessing it through random access") {
                            std::vector<std::size_t> mapping {3, 4};
                            auto ix = 0;
                            for(auto it = v.begin(); it != v.end(); ++it, ++ix) {
                                REQUIRE(it->val() == v.at(ix).val());
                            }
                        }
                        WHEN("Accessing it through plus operator") {
                            REQUIRE((*(v.begin() + 0)).val() == 8);
                            REQUIRE((*(v.cbegin() + 0)).val() == 8);
                            REQUIRE((*(v.begin() + 1)).val() == 3);
                            REQUIRE((*(v.cbegin() + 1)).val() == 3);
                        }
                        WHEN("Accessing it through minus operator") {
                            REQUIRE((v.end() - 1)->val() == 3);
                            REQUIRE((v.cend() - 1)->val() == 3);
                            REQUIRE((v.end() - 2)->val() == 8);
                            REQUIRE((v.cend() - 2)->val() == 8);
                        }
                        WHEN("Accessing it through stepwise increase") {
                            auto it = v.begin();
                            REQUIRE((it++)->val() == 8);
                            REQUIRE((it++)->val() == 3);
                            REQUIRE(it == v.end());

                            std::vector<int> vv;
                            assert(vv.begin() == vv.cend());
                        }
                        WHEN("Accessing it through stepwise increase const") {
                            auto it = v.begin();
                            REQUIRE((it++)->val() == 8);
                            REQUIRE((it++)->val() == 3);
                            REQUIRE(it == v.cend());
                        }
                        WHEN("Accessing it through stepwise decrease") {
                            auto itEnd = v.end();
                            REQUIRE((--itEnd)->val() == 3);
                            REQUIRE((--itEnd)->val() == 8);
                            REQUIRE(itEnd == v.begin());
                        }
                        WHEN("Accessing it through stepwise decrease const") {
                            auto itEnd = v.cend();
                            REQUIRE((--itEnd)->val() == 3);
                            REQUIRE((--itEnd)->val() == 8);
                            REQUIRE(itEnd == v.cbegin());
                        }
                    }

                    AND_WHEN("Removing element 4") {
                        v.erase(v.begin() + 1);
                        THEN("The active size is 1") {
                            REQUIRE(std::distance(v.begin(), v.end()) == 1);
                            REQUIRE(v.size() == 1);
                        }
                        AND_THEN("the active iterator skips the elements with numbers 7, 1, 5, 3") {
                            std::vector<std::size_t> mapping {3};
                            auto ix = 0;
                            for(auto it = v.begin(); it != v.end(); ++it, ++ix) {
                                REQUIRE(it->val() == v.at(ix).val());
                            }
                            REQUIRE(((v.begin() + 0))->val() == 8);
                            REQUIRE((v.end() - 1)->val() == 8);
                            REQUIRE((--v.end())->val() == 8);
                        }

                        AND_WHEN("Removing element 3") {
                            v.erase(v.begin());
                            THEN("The vector is empty") {
                                REQUIRE(v.empty());
                                REQUIRE(v.size() == 0);
                            }
                            AND_THEN("the active iterator skips all the elements") {
                                REQUIRE(v.begin() == v.end());
                            }
                        }

                        AND_WHEN("Adding a new element 100") {
                            v.push_back(A{100});
                            THEN("The active size is 2") {
                                REQUIRE(v.size() == 2);
                                REQUIRE(std::distance(v.begin(), v.end()) == 2);
                            }
                            AND_THEN("The vector contains {8, 100}") {
                                auto it8 = std::find_if(v.begin(), v.end(), [](auto a) { return a.val() == 8; });
                                REQUIRE(!it8->deactivated());
                                REQUIRE(it8 != v.end());
                                auto it100 = std::find_if(v.begin(), v.end(), [](auto a) { return a.val() == 100; });
                                REQUIRE(!it100->deactivated());
                                REQUIRE(it100 != v.end());
                            }
                        }
                    }
                }
            }
        }
    }
}
