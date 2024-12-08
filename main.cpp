#include <iostream>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <print>
#include <ctime>
#include <chrono>
#include <map>
#include <ranges>
#include <cassert>
#include <fstream>

// ONE UNIT = 1 YEAR
// YIELD = CONTINUOUS TIME YIELD
using TIME_POINT = double;
using AMOUNT = double;

using PaymentSchedule = std::vector<std::pair<TIME_POINT, AMOUNT>>;

enum class DayCountConvention
{
    _30_360,
    ACTUAL_ACTUAL
};

struct Bond
{
    Bond(double principal, double coupon, double frequency, double n_years, int start_date, DayCountConvention convention)
        : start_date(start_date), convention(convention)
    {
        assert(frequency > 0);
        assert(n_years >= 0);
        for (int i = 1; i <= n_years * frequency; ++i)
            schedule.emplace_back(i / frequency, coupon / frequency);
        schedule.back().second += principal;
    }

    [[nodiscard]] auto get_dirty_price_at_point(double point, double cts_yield) const noexcept
    {
        if (point >= schedule.back().first)
        {
            std::cerr << "Can't find bond price when it is already matured!";
            std::terminate();
        }

        const auto new_schedule = get_schedule_on_point(point);
        return std::accumulate(new_schedule.begin(), new_schedule.end(), 0.0, [cts_yield](double sum, auto const& p) {
            return sum + p.second * std::exp(-cts_yield * p.first);
        });
    }

    [[nodiscard]] auto get_clean_price_at_point(double point, double cts_yield) const noexcept
    {
        const auto dirty_price = get_dirty_price_at_point(point, cts_yield);
        const auto last_payment_point = std::find_if(schedule.rbegin(), schedule.rend(), [&](auto const& p) { return p.first <= point; });
        const auto start_point = (last_payment_point == schedule.rend()) ? 0 : last_payment_point->first;
        const auto interest = schedule[0].second * (point - start_point) / (schedule[1].first - schedule[0].first);
        return dirty_price - interest;
    }

    [[nodiscard]] auto get_dirty_price(int date, double cts_yield) const noexcept
    {
        return get_dirty_price_at_point(get_time_diff(start_date, date), cts_yield);
    }

    [[nodiscard]] auto get_clean_price(int date, double cts_yield) const noexcept
    {
        return get_clean_price_at_point(get_time_diff(start_date, date), cts_yield);
    }

    [[nodiscard]] auto get_cts_yield(int date, AMOUNT price) const noexcept
    {
        const auto payment_sum = std::accumulate(schedule.begin(), schedule.end(), 0.0, [](double sum, auto const& p) { return sum + p.second; });
        if (payment_sum < price)
        {
            std::cerr << "Price of bond can't be more than sum of payments. This implies negative yield which is not reasonable" << std::endl;
            std::terminate();
        }

        double l = 0;
        double r = 1;
        while (std::abs(l - r) > 1e-6)
        {
            const double m = std::midpoint(l, r);
            const double fetched_price = get_clean_price(date, m);
            if (fetched_price >= price)
                l = m;
            else
                r = m;
        }

        return l;
    }

private:
    /**
     * Returns the schedule after removing all payments before the given date. Note that the payment on passed point is
     * not included.
     * @param point The time point where units are measured in years, w.r.t. start date of bond.
     * @return The payment schedule, with timings reset to 0 w.r.t. passed point.
     */
    [[nodiscard]] PaymentSchedule get_schedule_on_point(double point) const noexcept
    {
        PaymentSchedule new_schedule;
        std::copy_if(schedule.begin(), schedule.end(), std::back_inserter(new_schedule), [point](auto const& p) { return p.first > point; });
        std::ranges::for_each(new_schedule, [point](auto& p) { p.first -= point; });
        return new_schedule;
    }

    /**
     * Finds the date difference in years w.r.t. day convention.
     * @param d1 First date.
     * @param d2 Second date.
     * @return Result is positive iff d1 <= d2.
     */
    [[nodiscard]] auto get_time_diff(int d1, int d2) const noexcept -> double
    {
        const int sign = d1 > d2 ? -1 : 1;
        if (sign == -1)
            std::swap(d1, d2);

        std::tm t1{}, t2{};
        t1.tm_year = (d1 / 10000) - 1900;
        t1.tm_mon = (d1 % 10000) / 100 - 1;
        t1.tm_mday = d1 % 100;

        t2.tm_year = (d2 / 10000) - 1900;
        t2.tm_mon = (d2 % 10000) / 100 - 1;
        t2.tm_mday = d2 % 100;

        if (convention == DayCountConvention::_30_360)
        {
            if (t2.tm_mday == 31 and (t1.tm_mday >= 30))
                t2.tm_mday = 30;
            if (t1.tm_mday == 31)
                t1.tm_mday = 30;
            return sign * ((t2.tm_year - t1.tm_year) * 360 + (t2.tm_mon - t1.tm_mon) * 30 + (t2.tm_mday - t1.tm_mday)) / 360.0;
        }

        if (convention == DayCountConvention::ACTUAL_ACTUAL)
        {
            const auto t1_seconds = std::mktime(&t1);
            const auto t2_seconds = std::mktime(&t2);
            return sign * std::abs(std::difftime(t2_seconds, t1_seconds) / (60 * 60 * 24 * 365.0));
        }

        std::cerr << "Unknown day count convention" << std::endl;
        std::terminate();
    }

    PaymentSchedule schedule;
    int start_date;
    DayCountConvention convention;
};

auto yearly_rate_to_cts(double rate, int frequency) -> double
{
    if (frequency == 0)
        return 0;

    return frequency * std::log(1 + rate / frequency);
}

auto cts_rate_to_yearly(double cts_rate, int frequency) -> double
{
    if (frequency == 0)
        return 0;

    return frequency * (std::exp(cts_rate / frequency) - 1);
}

int main()
{
    std::cout << std::setprecision(8) << std::fixed;
    const std::map<int, Bond> bonds = {
            {105759698, Bond(100, 3.6, 2, 15, 20190430, DayCountConvention::_30_360)},
            {165237838, Bond(100, 4, 2, 10, 20240208, DayCountConvention::ACTUAL_ACTUAL)},
            {168241282, Bond(100, 6.858, 2, 30, 20240501, DayCountConvention::_30_360)},
            {168302533, Bond(100, 4.625, 2, 30, 20240515, DayCountConvention::ACTUAL_ACTUAL)},
            {107424376, Bond(100, 3.75, 2, 30.5, 20190729, DayCountConvention::ACTUAL_ACTUAL)},
            {0, Bond(100, 2.4, 2, 10, 20160808, DayCountConvention::_30_360)}
    };

    const auto date = 20241101;
    const auto bond = bonds.at(168241282);
    const auto ref_bond = bonds.at(168302533);
    const auto bond_price = 107.056;
    const auto ref_price = 105.515625;

    auto bond_yield = cts_rate_to_yearly(bond.get_cts_yield(date, bond_price), 2);
    auto ref_yield = cts_rate_to_yearly(ref_bond.get_cts_yield(date, ref_price), 2);
    std::cout << "Bond Yield: " << bond_yield * 100 << "%" << std::endl;
    std::cout << "Ref Yield: " << ref_yield * 100 << "%" << std::endl;
    std::cout << "Spread: " << (bond_yield - ref_yield) * 10000 << "bps" << std::endl;

    return 0;
}
