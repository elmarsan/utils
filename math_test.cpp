#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <vector>
#include <iostream>

#include "math.h"

TEST_SUITE("Transformations")
{
    TEST_CASE("Uniform scaling: Magnification")
    {
        const auto s = scale(mat4{1.0f}, vec3{1.5});
        REQUIRE(s[0][0] == 1.5f);
        REQUIRE(s[1][1] == 1.5f);
        REQUIRE(s[2][2] == 1.5f);
        REQUIRE(s[3][3] == 1.0f);

        // Col0
        REQUIRE(s[0][1] == 0);
        REQUIRE(s[0][2] == 0);
        REQUIRE(s[0][3] == 0);
        // Col2
        REQUIRE(s[1][0] == 0);
        REQUIRE(s[1][2] == 0);
        REQUIRE(s[1][3] == 0);
        // Col2
        REQUIRE(s[2][0] == 0);
        REQUIRE(s[2][1] == 0);
        REQUIRE(s[2][3] == 0);
    }

    TEST_CASE("Uniform scaling: Reduction")
    {
        const auto s = scale(mat4{1.0f}, vec3{0.5});
        REQUIRE(s[0][0] == 0.5f);
        REQUIRE(s[1][1] == 0.5f);
        REQUIRE(s[2][2] == 0.5f);
        REQUIRE(s[3][3] == 1.0f);

        // Col0
        REQUIRE(s[0][1] == 0);
        REQUIRE(s[0][2] == 0);
        REQUIRE(s[0][3] == 0);
        // Col1
        REQUIRE(s[1][0] == 0);
        REQUIRE(s[1][2] == 0);
        REQUIRE(s[1][3] == 0);
        // Col2
        REQUIRE(s[2][0] == 0);
        REQUIRE(s[2][1] == 0);
        REQUIRE(s[2][3] == 0);
    }

    TEST_SUITE("Rotation")
    {
        auto testRotation = [](const mat4& m, const float a, const vec3& v) -> void
        {
            const auto c = cos(a);
            const auto s = sin(a);
            const auto r = v.normalize();
            const auto x = r.x;
            const auto y = r.y;
            const auto z = r.z;

            mat4 rot;
            rot[0][0] = ((1 - c) * pow(x, 2) + c);
            rot[0][1] = ((1 - c) * x * y + s * z);
            rot[0][2] = ((1 - c) * x * z - s * y);
            rot[0][3] = 0;

            rot[1][0] = ((1 - c) * x * y - s * z);
            rot[1][1] = ((1 - c) * pow(y, 2) + c);
            rot[1][2] = ((1 - c) * y * z + s * x);
            rot[1][3] = 0;

            rot[2][0] = ((1 - c) * x * z + s * y);
            rot[2][1] = ((1 - c) * y * z - s * x);
            rot[2][2] = ((1 - c) * pow(z, 2) + c);
            rot[2][3] = 0;

            rot[3][0] = 0;
            rot[3][1] = 0;
            rot[3][2] = 0;
            rot[3][3] = 0;

            const auto res = rotate(m, a, v);
            REQUIRE(res[0] == (rot[0][0] * m[0] + rot[0][1] * m[1] + rot[0][2] * m[2]));
            REQUIRE(res[1] == (rot[1][0] * m[0] + rot[1][1] * m[1] + rot[1][2] * m[2]));
            REQUIRE(res[2] == (rot[2][0] * m[0] + rot[2][1] * m[1] + rot[2][2] * m[2]));
            REQUIRE(res[3] == m[3]);
        };

        TEST_CASE("Identity matrix")
        {
            std::vector<vec3> rots{
                {1.0f, 0, 0},  // X-axis
                {0, 1.0f, 0},  // Y-axis
                {0, 0, 1.0f},  // Z-axis
            };

            // Random rotations
            for (int i = 0; i < 10; i++)
            {
                rots.push_back(vec3{
                    randomFloat(-1.0f, 2.0f),
                    randomFloat(-1.0f, 2.0f),
                    randomFloat(-1.0f, 2.0f),
                });
            }

            for (const auto rot : rots)
            {
                const auto a = randomFloat(0, 360.0f);
                testRotation(mat4{1.0f}, a, rot);
            }
        }

        TEST_CASE("Random matrix")
        {
            for (int i = 0; i < 10; i++)
            {
                const mat4 m{
                    // col0
                    randomFloat(-1.0f, 2.0f),
                    randomFloat(-1.0f, 2.0f),
                    randomFloat(-1.0f, 2.0f),
                    randomFloat(-1.0f, 2.0f),
                    // col1
                    randomFloat(-1.0f, 2.0f),
                    randomFloat(-1.0f, 2.0f),
                    randomFloat(-1.0f, 2.0f),
                    randomFloat(-1.0f, 2.0f),
                    // col2
                    randomFloat(-1.0f, 2.0f),
                    randomFloat(-1.0f, 2.0f),
                    randomFloat(-1.0f, 2.0f),
                    randomFloat(-1.0f, 2.0f),
                    // col3
                    randomFloat(-1.0f, 2.0f),
                    randomFloat(-1.0f, 2.0f),
                    randomFloat(-1.0f, 2.0f),
                    randomFloat(-1.0f, 2.0f),
                };

                const vec3 v{
                    randomFloat(-1.0f, 2.0f),
                    randomFloat(-1.0f, 2.0f),
                    randomFloat(-1.0f, 2.0f),
                };

                const auto a = radians(randomFloat(0, 360));
                testRotation(m, a, v);
            }
        }
    }
}

TEST_SUITE("Mat2")
{
    TEST_CASE("Multiplication")
    {
        const mat2 m1{
            3, 5,   // Col 0
            2, -3,  // Col 1
        };
        const mat2 m2{
            7, 1,   // Col 0
            1, -3,  // Col 1
        };
        const mat2 expected{
            23, 32,  // Col 0
            -3, 14   // Col 1
        };
        REQUIRE(m1 * m2 == expected);
    }

    TEST_CASE("Determinant")
    {
        std::vector<std::pair<mat2, float>> tests;
        // TEST 1
        tests.push_back({
            mat2{
                3, 5,   // Col 0
                -4, 7,  // Col 1
            },
            41,
        });
        // TEST 2
        tests.push_back({
            mat2{
                6, 7,  // Col 0
                8, 9,  // Col 1
            },
            -2,
        });
        // TEST 3
        tests.push_back({
            mat2{
                4, 7,  // Col 0
                7, 9,  // Col 1
            },
            -13,
        });
        // TEST 4
        tests.push_back({
            mat2{
                4, 6,  // Col 0
                7, 8,  // Col 1
            },
            -10,
        });
        // TEST 5
        tests.push_back({
            mat2{
                2, 3,  // Col 0
                8, 9,  // Col 1
            },
            -6,
        });
        // TEST 6
        tests.push_back({
            mat2{
                1, 3,  // Col 0
                7, 9,  // Col 1
            },
            -12,
        });

        for (int i = 0; i < tests.size(); i++)
        {
            const auto errMsg = "Failed test: " + std::to_string(i);
            REQUIRE_MESSAGE(tests[i].first.determinant() == tests[i].second, errMsg);
        }
    }

    TEST_CASE("Transpose")
    {
        const mat2 m{
            4, -3,  // Col 0
            8, 7,   // Col 1
        };
        const mat2 expected{
            4, 8,   // Col 0
            -3, 7,  // Col 1
        };
        REQUIRE(m.transpose() == expected);
    }

    TEST_CASE("Minor")
    {
        const mat2 m{
            1, 2,  // Col 0
            3, 4,  // Col 1
        };
        const mat2 expected{
            4, 3,  // Col 0
            2, 1,  // Col 1
        };
        REQUIRE(m.minor() == expected);
    }

    TEST_CASE("Cofactor")
    {
        const mat2 m{
            1, 2,  // Col 0
            3, 4,  // Col 1
        };
        const mat2 expected{
            4, -3,  // Col 0
            -2, 1,  // Col 1
        };
        REQUIRE(m.cofactor() == expected);
    }

    TEST_CASE("Adjugate")
    {
        const mat2 m{
            3, -2,  // Col 0
            7, -3,  // Col 1
        };
        const mat2 expected{
            -3, 2,  // Col 0
            -7, 3,  // Col 1
        };
        REQUIRE(m.adjugate() == expected);
    }

    TEST_CASE("Inverse")
    {
        const auto d = 41;
        const mat2 m{
            3, -7,  // Col 0
            5, 2,   // Col 1
        };
        const mat2 adjugate{
            2, 7,   // Col 0
            -5, 3,  // Col 1
        };
        REQUIRE(m.adjugate() == adjugate);
        REQUIRE(m.determinant() == d);

        const mat2 expected{
            adjugate[0][0] / d, adjugate[0][1] / d,  // Col 0
            adjugate[1][0] / d, adjugate[1][1] / d,  // Col 1
        };
        REQUIRE(m.inverse() == expected);

        const auto converse = m.inverse() * m;
        mat2 identity{1.0f};

        REQUIRE(converse[0][0] == doctest::Approx(identity[0][0]));
        REQUIRE(converse[0][1] == doctest::Approx(identity[0][1]));
        REQUIRE(converse[1][0] == doctest::Approx(identity[1][0]));
        REQUIRE(converse[1][1] == doctest::Approx(identity[1][1]));
    }
}

TEST_SUITE("Mat3")
{
    TEST_CASE("Multiplication")
    {
        const mat3 a{
            1,  3, -4,  // Col 0
            2,  2, 0,   // Col 1
            -1, 0, 2,   // Col 2
        };

        const mat3 b{
            3, 0, -2,  // Col 0
            4, 1, 0,   // Col 1
            2, 0, 1,   // Col 2
        };

        const mat3 expected{
            5, 9,  -16,  // Col 0
            6, 14, -16,  // Col 1
            1, 6,  -6,   // Col 2
        };

        REQUIRE(expected == a * b);
    }

    TEST_CASE("Determinant")
    {
        const mat3 m{
            5,  4,  1,   // Col 0
            7,  -3, 7,   // Col 1
            -8, 6,  -9,  // Col 2
        };
        REQUIRE(m.determinant() == -29.0f);

        const mat3 m2{
            2,  5, -8,  // Col 0
            4,  7, 1,   // Col 1
            -3, 6, 9,   // Col 2
        };
        REQUIRE(m2.determinant() == -441.0f);
    }

    TEST_CASE("Transpose")
    {
        const mat3 m{
            3,  1, 2,  // Col 0
            -5, 1, 3,  // Col 1
            -9, 2, 6,  // Col 2
        };
        const mat3 expected{
            3, -5, -9,  // Col 0
            1, 1,  2,   // Col 1
            2, 3,  6,   // Col 2
        };
        REQUIRE(m.transpose() == expected);
    }

    TEST_CASE("Minor")
    {
        std::vector<std::pair<mat3, mat3>> tests;
        // TEST 1
        tests.push_back({mat3{
                             1, 4, 7,  // Col 0
                             2, 6, 8,  // Col 1
                             3, 7, 9,  // Col 2
                         },
                         mat3{
                             -2, -6, -4,    // Col 0
                             -13, -12, -5,  // Col 1
                             -10, -6, -2,   // Col 2
                         }});
        // TEST 2
        tests.push_back({mat3{
                             3, 2, 7,    // Col 0
                             1, -3, -3,  // Col 1
                             2, 2, 0,    // Col 2
                         },
                         mat3{
                             6, 6, 8,       // Col 0
                             -14, -14, 2,   // Col 1
                             15, -16, -11,  // Col 2
                         }});

        for (int i = 0; i < tests.size(); i++)
        {
            const auto errMsg = "Failed test: " + std::to_string(i);
            REQUIRE_MESSAGE(tests[i].first.minor() == tests[i].second, errMsg);
        }
    }

    TEST_CASE("Cofactor")
    {
        const mat3 m{
            1, 4, 7,  // Col 0
            2, 6, 8,  // Col 1
            3, 7, 9,  // Col 2
        };
        const mat3 expected{
            -2,  6,   -4,  // Col 0
            13,  -12, 5,   // Col 1
            -10, 6,   -2,  // Col 2
        };
        REQUIRE(m.cofactor() == expected);
    }

    TEST_CASE("Adjugate")
    {
        const mat3 m{
            3,  7,  -1,  // Col 0
            -2, -3, 2,   // Col 1
            6,  8,  2,   // Col 2
        };
        const mat3 expected{
            -22, -22, 11,  // Col 0
            16,  12,  -4,  // Col 1
            2,   18,  5,   // Col 2
        };
        REQUIRE(m.adjugate() == expected);
    }

    TEST_CASE("Adjugate")
    {
        const mat3 m{
            3,  7,  -1,  // Col 0
            -2, -3, 2,   // Col 1
            6,  8,  2,   // Col 2
        };
        const mat3 expected{
            -22, -22, 11,  // Col 0
            16,  12,  -4,  // Col 1
            2,   18,  5,   // Col 2
        };
        REQUIRE(m.adjugate() == expected);
    }

    TEST_CASE("Inverse")
    {
        const mat3 m{
            4,  4,  1,   // Col 0
            5,  3,  -5,  // Col 1
            -5, -4, 3,   // Col 2
        };
        const mat3 expected{
            1.22f,  1.77f,  2.55f,   // Col 0
            -1.11f, -1.88f, -2.77f,  // Col 1
            0.55f,  0.44f,  0.88f,   // Col 2
        };

        const auto inverse = m.inverse();

        const float epsilon = 0.01f;
        for (int c = 0; c < 3; c++)
        {
            for (int r = 0; r < 3; r++)
            {
                const auto roundedVal = std::round(inverse[c][r] * 100.0f) / 100.0f;
                const auto roundedExpected = std::round(expected[c][r] * 100.0f) / 100.0f;
                REQUIRE(std::abs(roundedVal - roundedExpected) < epsilon);
            }
        }
    }
}

TEST_SUITE("Mat4")
{
    TEST_CASE("Determinant")
    {
        const mat4 m{
            3,  1,  3,  -6,  // Col 0
            -7, -4, 2,  6,   // Col 1
            8,  -3, 5,  2,   // Col 2
            1,  8,  -5, 1,   // Col 3
        };
        REQUIRE(m.determinant() == 2322);
    }

    /* TEST_CASE("Transpose") */
    /* { */
    /*     const mat4 m{ */
    /*         2,  4, 3, 3,  // Col 0 */
    /*         3,  2, 9, 7,  // Col 1 */
    /*         10, 2, 2, 7,  // Col 2 */
    /*         8,  1, 3, 1,  // Col 3 */
    /*     }; */
    /*     const mat4 expected{ */
    /*         2, 3, 10, 8,  // Col 0 */
    /*         4, 2, 2,  1,  // Col 1 */
    /*         3, 9, 2,  3,  // Col 2 */
    /*         3, 7, 7,  1,  // Col 3 */
    /*     }; */
    /*     REQUIRE(expected == m.transpose()); */
    /*     REQUIRE(m.transpose().transpose() == m); */
    /* } */

    /* TEST_CASE("Adjugate") */
    /* { */
    /*     const mat4 m{ */
    /*         -2, 4, -3, 2,  // Col 0 */
    /*         5,  1, 5,  2,  // Col 1 */
    /*         1,  0, 5,  3,  // Col 2 */
    /*         5,  3, 1,  3,  // Col 3 */
    /*     }; */
    /*     const mat4 expected{ */
    /*         27,  -72,  -27, 36,   // Col 0 */
    /*         9,   -108, -72, 117,  // Col 1 */
    /*         15,  51,   6,   -78,  // Col 2 */
    /*         -39, 69,   60,  -87,  // Col 3 */
    /*     }; */
    /*     REQUIRE(expected == m.adjugate()); */

    /*     const mat4 m2{ */
    /*         3,  -8, 3, 7,   // Col 0 */
    /*         -2, 7,  4, -2,  // Col 1 */
    /*         9,  7,  8, 2,   // Col 2 */
    /*         2,  2,  0, 6,   // Col 3 */
    /*     }; */

    /*     const mat4 expected2{ */
    /*         216,  264,  -434, -160,  // Col 0 */
    /*         610,  -118, -542, -164,  // Col 1 */
    /*         -386, -40,  -52,  142,   // Col 2 */
    /*         80,   -334, 343,  -563,  // Col 3 */
    /*     }; */
    /*     REQUIRE(expected == m.adjugate()); */
    /* } */

    /* TEST_CASE("Inverse") */
    /* { */
    /*     mat4 m{ */
    /*         3,  -8, 3, 7,   // Col 0 */
    /*         -2, 7,  4, -2,  // Col 1 */
    /*         9,  7,  8, 2,   // Col 2 */
    /*         2,  2,  0, 6,   // Col 3 */
    /*     }; */

    /*     mat4 expected{ */
    /*         -0.0555f, -0.1569f, 0.0993f,    -0.0205f,  // Col 0 */
    /*         -0.0679f, 0.0303f,  0.0102f,    0.0859f,   // Col 1 */
    /*         0.1116f,  0.1394f,  0.0133814f, -0.0882f,  // Col 2 */
    /*         0.0411f,  0.04220f, -0.03654f,  0.1448f,   // Col 3 */
    /*     }; */

    /*     const auto i = m.inverse(); */

    /*     for (int c = 0; c < 4; c++) */
    /*     { */
    /*         for (int r = 0; r < 4; r++) */
    /*         { */
    /*             REQUIRE(expected[c][r] == doctest::Approx(i[c][r]).epsilon(.05)); */
    /*         } */
    /*     } */
    /* } */
}
