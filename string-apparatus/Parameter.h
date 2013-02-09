// simple lerper class
class Parameter {
  Parameter() : targetValue(0.0f),
		previousValue ( 0.0f )
  {}
  inline double getValueAt(double u) {
    return lastValue = u * targetValue + (1.0f-u) * previousValue;
  }
  inline void setTargetValue( double target ) {
    targetValue = target;
    previousValue = lastValue;
  }
  inline void setCurrentValue ( double v ) {
    previousValue = lastValue = v;
  }
  double targetValue;
  double previousValue;
  double lastValue;
};
