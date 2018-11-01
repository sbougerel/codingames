#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <cstdio>
#include <thread>
#include <array>
#include <cmath>

using namespace std;

const int GAME_RAD = 6000;
const int WATERTOWN_RAD = 3000;
const int SKILL_RANGE = 2000;
const int SKILL_RAD = 1000;
const int SKILL_COST = 30;
const string WAIT = "WAIT";
const double REAPER_FRICTION = 0.2;
const double DESTROYER_FRICTION = 0.3;
const double TANKER_FRICTION = 0.4;

enum UnitType : int {
    Reaper,
    Destroyer,
    Tanker,
    Wreck,
    Tarpool,
    Oilpool,
    Ignore,
};

enum PlayerId : int {
    Me,
    Them,
};

template<typename Tp>
struct Vec2 : array<Tp, 2>
{
    typedef array<Tp, 2> Base;

    Vec2() {}
    Vec2(Tp a, Tp b)
    { Base::operator[](0) = a; Base::operator[](1) = b; }

    Tp x() const { return Base::operator[](0); }
    Tp& x() { return Base::operator[](0); }

    Tp y() const { return Base::operator[](1); }
    Tp& y() { return Base::operator[](1); }
};

template<typename Tp>
inline Vec2<Tp> operator+ (const Vec2<Tp>& a, const Vec2<Tp>& b)
{ return {a.x() + b.x(), a.y() + b.y()}; }

template<typename Tp>
inline Vec2<Tp> operator- (const Vec2<Tp>& a, const Vec2<Tp>& b)
{ return {a.x() - b.x(), a.y() - b.y()}; }

template<typename Tp>
inline Vec2<Tp> operator- (const Vec2<Tp>& a)
{ return {-a.x(), -a.y()}; }

template<typename Tp>
inline Vec2<Tp> operator* (const Vec2<Tp>& a, Tp factor)
{ return {a.x() * factor, a.y() * factor}; }

template<typename Tp>
inline Vec2<Tp> operator/ (const Vec2<Tp>& a, Tp factor)
{ return {a.x() / factor, a.y() / factor}; }

template<typename Tp>
constexpr Tp sq(Tp x) { return x * x; }

template<typename Tp>
constexpr Tp dist(const Vec2<Tp>& a, const Vec2<Tp>& b)
{ return sq(a.x() - b.x()) + sq(a.y() - b.y()); } // won't overflow, so we're good

template<typename Tp>
constexpr Tp magsq(const Vec2<Tp>& a)
{ return sq(a.x()) + sq(a.y()); }

inline double mag(const Vec2<double>& a)
{ return sqrt(sq(a.x()) + sq(a.y())); }

inline Vec2<double> normalize(const Vec2<double>& a, double norm = 1.)
{
    double ratio = norm / mag(a);
    return {a.x() * ratio, a.y() * ratio};
}

// sqrt can't be computed on integer, and there's a fast convergent approx with 2 iterations only
inline int mag(const Vec2<int>& a)
{
    // the squared value
    int S = sq(a.x()) + sq(a.y());
    // start at Manhattan magnitude
    int x = abs(a.x()) + abs(a.y());
    // do 2 increments of newton's method, add +1 to avoid div by 0
    x = (sq(x) + S) / (2 * x + 1);
    return (sq(x) + S) / (2 * x + 1);
}

inline Vec2<int> normalize(const Vec2<int>& a, int norm)
{
    int x = mag(a);
    // avoid overflows/underflows, and div by 0
    return {(a.x() * norm) / (x + 1), (a.y() * norm) / (x + 1)};
}

inline bool collide(const Vec2<int> a, const Vec2<int> b, int rad1, int rad2)
{ return ((sq(rad1) + sq(rad2)) > dist(a, b)); }


struct Unit {
    Unit() : id(0), type(UnitType::Ignore), player(PlayerId::Me),
        mass(.0), radius(0), pos{0, 0}, speed{0, 0}, accel{0, 0},
        extra(-1), extra2(-1)
        { }

    int id;
    int type;
    int player;
    double mass;
    int radius;
    Vec2<int> pos;
    Vec2<int> speed;
    Vec2<int> accel;
    int extra;
    int extra2;

    int x() const { return pos.x(); }
    int& x() { return pos.x(); }

    int y() const { return pos.y(); }
    int& y() { return pos.y(); }

    int vx() const { return speed.x(); }
    int& vx() { return speed.x(); }

    int vy() const { return speed.y(); }
    int& vy() { return speed.y(); }

    int ax() const { return accel.x(); }
    int& ax() { return accel.x(); }

    int ay() const { return accel.y(); }
    int& ay() { return accel.y(); }

    int waterCapacity() const { return extra2; } // aliases
    int waterQty() const { return extra; } // aliases
    int duration() const { return extra; } // aliases
};

typedef vector<Unit> Units;

struct Frame {
    Frame() : myScore(0), theirScore(0), myRage(0), theirRage(0),
        myReap(nullptr), theirReap(nullptr), myDestroy(nullptr), theirDestroy(nullptr),
        units() { }

    int myScore;
    int theirScore;
    int myRage;
    int theirRage;
    Unit* myReap;
    Unit* theirReap;
    Unit* myDestroy;
    Unit* theirDestroy;
    Units units;
};

inline istream& operator >> (istream& in, Frame& fr) {
    in >> fr.myScore; in.ignore();
    in >> fr.theirScore; in.ignore();
    in >> fr.myRage; in.ignore();
    in >> fr.theirRage; in.ignore();
    size_t unitCount;
    in >> unitCount; in.ignore();

    // Reuse memory
    fr.units.reserve(unitCount);
    for (size_t i = 0; i < unitCount; i++) {
        if ( fr.units.size() == i ) { fr.units.push_back(Unit()); };
        Unit &u = fr.units[i];
        in >> u.id >> u.type >> u.player >> u.mass >> u.radius >> u.x() >> u.y() >> u.vx() >> u.vy() >> u.extra >> u.extra2;
        in.ignore();

        // organise stuffs
        if (u.player == PlayerId::Me) {
            switch (u.type) {
            case UnitType::Destroyer:
                fr.myDestroy = &u; break;
            case UnitType::Reaper:
                fr.myReap = &u; break;
            default:
                cerr << "Can't find my units!" << endl;
            }
        }
        if (u.player == PlayerId::Them) {
            switch (u.type) {
            case UnitType::Destroyer:
                fr.theirDestroy = &u; break;
            case UnitType::Reaper:
                fr.theirReap = &u; break;
            default:
                cerr << "Can't find my units!" << endl;
            }
        }
    }

    // Ignore units in extra memory
    for (size_t i = unitCount; i < fr.units.size(); ++i) {
        fr.units[i].type = UnitType::Ignore;
    }
    return in;
}

inline void derive(Frame& curr, const Frame& prev)
{
    curr.myDestroy->accel = curr.myDestroy->speed - prev.myDestroy->speed;
    curr.theirDestroy->accel = curr.theirDestroy->speed - prev.theirDestroy->speed;
    curr.myReap->accel = curr.myReap->speed - prev.myReap->speed;
    curr.theirReap->accel = curr.theirReap->speed - prev.theirReap->speed;
}

inline const Unit* find_collision(const Unit& a, const Units& us, UnitType ut, int rad)
{
    for (const Unit& u : us) {
        if (u.type == ut && collide(a.pos, u.pos, u.radius, a.radius))
            return &u;
    }
    return nullptr;
}

// Commands
inline
void Grenade(const Vec2<int>& pos, string msg = string()) {
    cout << "GRENADE " << pos.x() << " " << pos.y() << " " << msg << endl;
}

inline
void Tar(const Vec2<int>& pos, string msg = string()) {
    cout << "TAR " << pos.x() << " " << pos.y() << " " << msg << endl;
}

inline
void Oil(const Vec2<int>& pos, string msg = {}) {
    cout << "OIL " << pos.x() << " " << pos.y() << " " << msg << endl;
}

inline void Ram(const Unit& ram, const Unit& target, string msg = {}) {
    // When ramming, ignore our speed, consider target speed over distance
    int distsq = dist(ram.pos, target.pos);
    int mssq = magsq(ram.speed);
    int ratio = distsq / mssq; // ratio of extrapolation
    auto proj = target.pos + (target.speed + target.accel / 2) * ratio;
    cout << proj.x() << " " << proj.y() << " 300 " << msg << endl;
}

inline void Reach(const Unit& unit, const Vec2<int>& pos, string msg = {}) {
    // Adjust "toward" direction based on speed.
    auto toward = pos - unit.speed;
    // Adjust accelertion based on distance
    auto toward_speed_delta = (toward - unit.pos) - unit.speed;
    int throttle = magsq(toward_speed_delta) / sq(unit.mass);
    cout << toward.x() << " " << toward.y() << " " << throttle << " " << msg << endl;
}

// A Priority is a pointer to valuable unit, a heuristic score (min - max)
// Heuristic is = waterQty * (mydistance - theirdistance)
// the smaller the value the better for us.
struct Priority {
    const Unit* unit;
    int heuristic;
};

// Heuristic is = (mydistance / waterQty)
// the smaller the value the better for us.
inline int heuristic(const Unit& target, const Unit& myDestroyer) {
    return dist(target.pos, myDestroyer.pos) / target.waterQty();
}

inline bool operator< (const Priority& a, const Priority& b) {
    return (a.heuristic < b.heuristic);
}

typedef vector<Priority> Priorities;

struct ReaperTactics {
    // The reaper works more like a state machine, cycling through various actions:
    enum ReapState : int {
        Annoy,
        EvictReap,
        ComboReap,
        EvictDestroy,
        ComboDestroy
    };

    ReaperTactics(): state(ReapState::Annoy) { }

    // Go through the state machine
    void updateAction(const Frame& curr, const Frame& prev)
    {
        switch (state) {
            case ReapState::Annoy: {
                if (curr.myRage > SKILL_COST * 2) { // so we can combo!
                    // if ennemy is close to a wreck...
                    for (const Unit& u : curr.units) {
                        if (u.type == UnitType::Wreck
                            && dist(curr.theirDestroy->pos, u.pos) < sq(SKILL_RAD + u.radius)) {
                            state = ReapState::EvictDestroy;
                            break;
                        }
                    }
                    // if their reap has been bothering our destroy
                    if (state == ReapState::Annoy && prev.myDestroy != nullptr
                        && dist(curr.myDestroy->pos, curr.theirReap->pos) < sq(SKILL_RAD + curr.myDestroy->radius)
                        && dist(prev.myDestroy->pos, prev.theirReap->pos) < sq(SKILL_RAD + curr.myDestroy->radius)) {
                        state = ReapState::EvictReap;
                    }
                    else if (curr.myRage > SKILL_COST * 3) // really need to do something :)
                        state = ReapState::EvictDestroy;
                }
                Ram(*curr.myReap, *curr.theirDestroy, "Ram");
            }
            break;
            case ReapState::EvictReap: {
                // Grenade my own destroy to evict their pesky reap
                if (dist(curr.myReap->pos, curr.myDestroy->pos) < sq(SKILL_RANGE)) {
                    state = ReapState::ComboReap;
                    Grenade(curr.myDestroy->pos, " Push Reap");
                }
                else {
                    Ram(*curr.myReap, *curr.myDestroy, "Push Reap");
                }
            }
            break;
            case ReapState::ComboReap: {
                // If their reap is far enough from my destroyer, put a tar pool right under it.
                if (dist(curr.theirReap->pos, curr.myDestroy->pos) > sq(SKILL_RAD + curr.myDestroy->radius)
                    && dist(curr.myReap->pos, curr.theirReap->pos) < sq(SKILL_RANGE)) {
                    state = ReapState::Annoy;
                    Tar(curr.theirReap->pos, "Trap Reap");
                }
                else {
                    state = ReapState::Annoy;
                    Ram(*curr.myReap, *curr.theirReap, "Ram");
                }
            }
            break;
            case ReapState::EvictDestroy: {
                // Try to push away destroy from center:
                auto pos = curr.theirDestroy->pos + normalize(-curr.theirDestroy->pos, SKILL_RAD);
                // Unless our own destroy is pushed
                if (dist(curr.myDestroy->pos, pos) > sq(SKILL_RAD + curr.myDestroy->radius)) {
                    if (dist(curr.myReap->pos, pos) > sq(SKILL_RANGE)) {
                        Reach(*curr.myReap, pos, "Push Destroy");
                    }
                    else {
                        state = ReapState::ComboDestroy;
                        Grenade(pos, "Push Destroy");
                    }
                }
                else {
                    // Center grenade on our destroy instead
                    if (dist(curr.myReap->pos, curr.myDestroy->pos) > sq(SKILL_RANGE)) {
                        Ram(*curr.myReap, *curr.myDestroy, "Push Destroy");
                    }
                    else {
                        state = ReapState::ComboDestroy;
                        Grenade(curr.myDestroy->pos, "Push Destroy");
                    }
                }
            }
            break;
            case ReapState::ComboDestroy: {
                // If their destroy is far enough from my destroyer, put a tar pool right under it.
                if (dist(curr.theirDestroy->pos, curr.myDestroy->pos) > sq(SKILL_RAD + curr.myDestroy->radius)
                    && dist(curr.myReap->pos, curr.theirDestroy->pos) < sq(SKILL_RANGE)) {
                    state = ReapState::Annoy;
                    Tar(curr.theirDestroy->pos, "Trap Destroy");
                }
                else {
                    state = ReapState::Annoy;
                    Ram(*curr.myReap, *curr.theirDestroy, "Trap Destroy");
                }
            }
        }
    }

    ReapState state;
};

struct DestroyerTactics {
    void updatePriorities(const Frame& curr, const Frame& prev)
    {
        priorities.resize(0); // remove elements but don't deallocate
        priorities.reserve(curr.units.size());
        // Establish the list of priorities:
        for (const Unit& u : curr.units) {
            switch (u.type) {
                case UnitType::Tanker: {
                    // if the tanker is approaching, good.
                    int mag_speed = mag(u.speed);
                    if (mag(u.speed + normalize(u.pos, mag_speed)) < mag_speed) {
                        // if it's within bounds, add it, otherwise don't care.
                        if (magsq(u.pos) < sq(GAME_RAD))
                            priorities.push_back({&u, heuristic(u, *curr.myDestroy)});
                    }
                    else {
                        // If still within game radius but closer to exit that we are to it
                        if (sq(u.radius) + magsq(u.pos) < sq(GAME_RAD)
                            && sq(GAME_RAD) - magsq(u.pos) > dist(u.pos, curr.myDestroy->pos))
                            priorities.push_back({&u, heuristic(u, *curr.myDestroy)});
                    }
                }
                break;
                case UnitType::Wreck: {
                    // if it is covered by an oil pool, ignore it.
                    priorities.push_back({&u, heuristic(u, *curr.myDestroy)});
                }
                default:
                break;
            }
        }
        sort(priorities.begin(), priorities.end(), [](const Priority& a, const Priority& b){ return (a < b); });
    }

    // Destroyer destroys highest priority target
    void updateAction(const Frame& curr, const Frame& prev)
    {
        updatePriorities(curr, prev);

        // if I'm in an wreck, collect it!
        const Unit* wreck = find_collision(*curr.myDestroy, curr.units, UnitType::Wreck, curr.myDestroy->radius);
        if (wreck != nullptr) {
            cout << wreck->x() << " " << wreck->y() << " 300 " << " take: " << wreck->id << endl;
        }
        else if (priorities.size() > 0) {
            const Unit* destroy_target = priorities[0].unit;
            cout << destroy_target->x() << " " << destroy_target->y() << " 300 " << " aim: " << destroy_target->id << endl;
        }
        else {
            // annoys enemy reaper otherwise
            const Unit* destroy_target = curr.theirDestroy;
            cout << destroy_target->x() << " " << destroy_target->y() << " 300" << endl;
        }
    }

    vector<Priority> priorities;
};

struct Game {
    Game() : frames(), curr(&frames[0]), prev(&frames[1]), rTactics(), dTactics() { }

    void updateFrames(istream& in) {
        swap(curr, prev);
        in >> *curr;
        if (prev->units.size() > 0) derive(*curr, *prev); // except first round
    }

    void run() {
        rTactics.updateAction(*curr, *prev);
        dTactics.updateAction(*curr, *prev);
    }

    array<Frame, 2> frames; // store frames
    Frame *curr; // shorthand for the current frame
    Frame *prev; // shorthand for the last frame
    ReaperTactics rTactics;
    DestroyerTactics dTactics;
};

/**
 * Auto-generated code below aims at helping you parse
 * the standard input according to the problem statement.
 **/
int main()
{
    Game game;

    // game loop
    while (1) {
        game.updateFrames(cin);
        game.run();
    }
}
