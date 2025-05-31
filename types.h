#ifndef TYPES_H
#define TYPES_H

struct ScoreTime {
    int score = 0;
    float time = 0;
    ScoreTime() = default;
    ScoreTime(int s, double t) : score(s), time(t) {}
};


#endif // TYPES_H
