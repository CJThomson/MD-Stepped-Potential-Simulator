#pragma once
struct eventTimes
{
  typedef enum {
    IP_OUT = 0,
    IP_IN = 1,
    SENTINAL = 2,
    NEIGHBOURCELL = 3,
    WALL = 4,
    NONE = 5
  } EventType;

eventTimes(double dt = HUGE_VAL, int p1 = -1, int p2 = -1, int coll = 0, EventType type = NONE):
  collisionTime(dt),
    particle1(p1),
    particle2(p2),
    p2coll(coll),
    _type(type)
  {}

  double collisionTime;
  int particle1;
  int particle2;
  int p2coll;
  EventType _type;

  bool operator<(const eventTimes& otherevent) const
  { return collisionTime < otherevent.collisionTime; }
  bool operator>(const eventTimes& otherevent) const
  { return collisionTime > otherevent.collisionTime; }
};
