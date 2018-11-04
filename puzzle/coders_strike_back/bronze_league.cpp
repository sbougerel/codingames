#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <array>
#include <cmath>

using namespace std;

template<typename Tp>
constexpr Tp sq(Tp x) { return x * x; }

template<typename Tp>
inline Tp cosi(Tp a, Tp unit) {
    double angle = (double)a;
    return static_cast<Tp>(cos(angle*M_PI/180.0) * (double)unit);
}

inline void thrust(int x, int y, int t) {
    if (t > 100) t = 100;
    if (t < 0) t = 0;
    cout << x << " " << y << " " << t << endl;
}

inline void boost(int x, int y) {
    cout << x << " " << y << " BOOST" << endl;
}


/**
 * Auto-generated code below aims at helping you parse
 * the standard input according to the problem statement.
 **/
int main()
{
    bool boost_used = false;

    // game loop
    while (1) {
        int x;
        int y;
        int nextCheckpointX; // x position of the next check point
        int nextCheckpointY; // y position of the next check point
        int nextCheckpointDist; // distance to the next checkpoint
        int nextCheckpointAngle; // angle between your pod orientation and the direction of the next checkpoint
        cin >> x >> y >> nextCheckpointX >> nextCheckpointY >> nextCheckpointDist >> nextCheckpointAngle; cin.ignore();
        int opponentX;
        int opponentY;
        cin >> opponentX >> opponentY; cin.ignore();

        if (nextCheckpointDist < 600) {
            if (abs(nextCheckpointAngle) < 70)
                thrust(nextCheckpointX, nextCheckpointY, cosi(nextCheckpointAngle, 100) - 20);
            else
                thrust(nextCheckpointX, nextCheckpointY, 0);
        }
        else if (nextCheckpointDist < 1000) {
            if (abs(nextCheckpointAngle) < 90)
                thrust(nextCheckpointX, nextCheckpointY, cosi(nextCheckpointAngle, 100));
            else
                thrust(nextCheckpointX, nextCheckpointY, 0);
        }
        else if (nextCheckpointDist < 3000) {
            if (abs(nextCheckpointAngle) < 100)
                thrust(nextCheckpointX, nextCheckpointY, cosi(nextCheckpointAngle, 100) + 30);
            else
                thrust(nextCheckpointX, nextCheckpointY, 0);
        }
        else {
            if (abs(nextCheckpointAngle) < 30 && !boost_used) {
                boost(nextCheckpointX, nextCheckpointY);
                boost_used = true;
            }
            else if (abs(nextCheckpointAngle) < 100)
                thrust(nextCheckpointX, nextCheckpointY, cosi(nextCheckpointAngle, 100) + 50);
            else
                thrust(nextCheckpointX, nextCheckpointY, 0);
        }
    }
}
