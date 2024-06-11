#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <vector>

#include "math.h"

TEST_SUITE("Transformations")
{
    TEST_CASE("Uniform scaling: Magnification")
    {
        auto s = scale(mat4{1.0f}, vec3{1.5});
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
        auto s = scale(mat4{1.0f}, vec3{0.5});
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
                mat4 m{
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

                vec3 v{
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
