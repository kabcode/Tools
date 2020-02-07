#ifndef HighPrecisionTimer_H
#define HighPrecisionTimer_H

#include <chrono>
#include <ostream>

using ClockType = std::conditional_t<
                    std::chrono::high_resolution_clock::is_steady,
                    std::chrono::high_resolution_clock,
                    std::chrono::steady_clock>;

enum class TimeUnits
{
    Nanoseconds,
    Microseconds,
    Milliseconds,
    Seconds
};

template<TimeUnits unit, bool Print = true>
class HighPrecisionTimer
{

public:
    HighPrecisionTimer(std::ostream& os = std::cout): m_os(os) 
    {
        
    };

    HighPrecisionTimer(long long* ptr_extern, std::ostream& os = std::cout) : m_os(os)
    {
        ptr_duration = ptr_extern;
    }

    inline ~HighPrecisionTimer()
    {
        ClockType::time_point Stop = ClockType::now();
        if (unit == TimeUnits::Nanoseconds)
            m_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop - m_Start).count();
        if (unit == TimeUnits::Microseconds)
            m_duration = std::chrono::duration_cast<std::chrono::microseconds>(Stop - m_Start).count();
        if (unit == TimeUnits::Milliseconds)
            m_duration = std::chrono::duration_cast<std::chrono::milliseconds>(Stop - m_Start).count();
        if (unit == TimeUnits::Seconds)
            m_duration = std::chrono::duration_cast<std::chrono::seconds>(Stop - m_Start).count();

        if constexpr (Print == true) {
            if constexpr (unit == TimeUnits::Nanoseconds) { m_timeunit = "ns"; }
            if constexpr (unit == TimeUnits::Microseconds) { m_timeunit = "us"; }
            if constexpr (unit == TimeUnits::Milliseconds) { m_timeunit = "ms"; }
            if constexpr (unit == TimeUnits::Seconds) { m_timeunit = "s"; }
             m_os << m_duration << " " << m_timeunit << std::endl;
        }

        if(ptr_duration != nullptr)
            *ptr_duration = m_duration;
    };

    static TimeUnits getTimeUnit() { return unit; }

private:
    long long m_duration = 0;
    std::string m_timeunit;
    std::ostream& m_os;
    long long* ptr_duration = nullptr;
    ClockType::time_point m_Start = ClockType::now();
};

#endif
