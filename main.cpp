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

enum class RateType
{
    YEARLY,
    CONTINUOUS
};

struct Rate
{
private:
    double cts_rate;

    Rate(double rate, RateType type) : cts_rate(type == RateType::YEARLY ? std::log(1 + rate) : rate) {
        if (rate > 1 || rate < 0)
        {
            std::cerr << "Rate should be between 0 and 1" << std::endl;
            std::terminate();
        }
    }

public:
    static Rate make_continuous_rate(double rate) noexcept
    {
        return {rate, RateType::CONTINUOUS};
    }

    static Rate make_yearly_rate(double rate, int frequency) noexcept
    {
        return {std::pow(1 + rate / frequency, frequency) - 1, RateType::YEARLY};
    }

    [[nodiscard]] auto get_continuous_rate() const noexcept
    {
        return cts_rate;
    }

    [[nodiscard]] auto get_yearly_rate(int frequency = 1) const noexcept
    {
        if (frequency == 0)
            return 0.0;

        return frequency * (std::exp(cts_rate / frequency) - 1);
    }

    [[nodiscard]] auto get_interest(double amount, double time, RateType rate_type) const noexcept
    {
        if (time < 0)
        {
            std::cerr << "Time should be non-negative" << std::endl;
            std::terminate();
        }

        if (rate_type == RateType::CONTINUOUS)
            return amount * (std::exp(cts_rate * time) - 1);

        if (int(time) != time and time > 1)
        {
            std::cerr << "Time should be an integer, or less than 1 for yearly rate. Passed: " << time << std::endl;
            std::terminate();
        }

        const auto yearly_rate = get_yearly_rate(1);

        if (time < 1)
            return amount * yearly_rate * time;

        return amount * std::pow(1 + yearly_rate, time) - amount;
    }

    [[nodiscard]] auto get_discount_factor(double time, RateType rate_type) const noexcept
    {
        return 1.0 / (1 + get_interest(1, time, rate_type));
    }

    [[nodiscard]] auto get_present_value(double amount, double time, RateType rate_type) const noexcept
    {
        return amount * get_discount_factor(time, rate_type);
    }

    [[nodiscard]] auto get_future_value(double amount, double time, RateType rateType) const noexcept
    {
        return amount / get_discount_factor(time, rateType);
    }

    [[nodiscard]] static bool test_rate_api() noexcept
    {
        if (std::abs(Rate(0.05, RateType::CONTINUOUS).get_continuous_rate() - 0.05) > 1e-6)
        {
            std::cerr << "Continuous rate fetch failed" << std::endl;
            return false;
        }

        const auto rate = Rate(0.05, RateType::YEARLY);
        if (std::abs(rate.get_yearly_rate(1) - 0.05) > 1e-6)
        {
            std::cerr << "Yearly rate fetch failed" << std::endl;
            return false;
        }

        if (std::abs(rate.get_interest(100, 1, RateType::CONTINUOUS) - 5) > 1e-6)
        {
            std::cerr << "Continuous interest calculation failed" << std::endl;
            return false;
        }

        if (std::abs(rate.get_interest(100, 1, RateType::YEARLY) - 5) > 1e-6)
        {
            std::cerr << "Yearly interest calculation failed" << std::endl;
            return false;
        }

        if (std::abs(rate.get_interest(100, 0.5, RateType::YEARLY) - 2.5) > 1e-6)
        {
            std::cerr << "Yearly interest calculation failed for non-integer time" << std::endl;
            return false;
        }

        if (std::abs(rate.get_interest(100, 0.5, RateType::CONTINUOUS) - 2.4695076) > 1e-6)
        {
            std::cerr << "Yearly interest calculation failed for non-integer time" << std::endl;
            return false;
        }

        if (std::abs(rate.get_discount_factor(1, RateType::YEARLY) - 0.95238095238) > 1e-6)
        {
            std::cerr << "Yearly discount factor calculation failed" << std::endl;
            return false;
        }

        if (std::abs(rate.get_discount_factor(2, RateType::YEARLY) - 0.90702947845) > 1e-6)
        {
            std::cerr << "Yearly discount factor calculation failed" << std::endl;
            return false;
        }

        if (std::abs(rate.get_discount_factor(2, RateType::CONTINUOUS) - 0.90702947845) > 1e-6)
        {
            std::cerr << "Yearly discount factor calculation failed" << std::endl;
            return false;
        }

        return true;
    }
};

enum class DayCountConvention
{
    _30_360,
    ACTUAL_ACTUAL
};

/**
* Finds the date difference in years w.r.t. day convention.
* @param d1 First date.
* @param d2 Second date.
* @return Result is positive iff d1 <= d2.
*/
[[nodiscard]] auto get_time_diff(DayCountConvention convention, int d1, int d2)
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

class Bond
{
public:
    double principal;
    double yearly_coupon;
    double n_years;
    int start_date;
    int frequency;
    DayCountConvention convention;

    [[nodiscard]] auto get_price(int date, Rate yield) const noexcept
    {
        if (int(n_years * frequency) != n_years * frequency)
        {
            std::cerr << "Bond should have integer number of payments" << std::endl;
            std::terminate();
        }

        const auto point = get_time_diff(convention, start_date, date);
        const auto interval_start = std::floor(point * frequency) / (double)frequency;
        const auto interval_end = interval_start + 1 / (double)frequency;
        if (point < 0 or point >= n_years)
        {
            std::cerr << "Can't find bond price when it is already matured!";
            std::terminate();
        }

        struct
        {
            double clean_price = 0;
            double dirty_price = 0;
        } price;

        price.dirty_price = 0.0;
        for (int i = (int)std::floor(point * frequency) + 1; i <= n_years * frequency; ++i)
            price.dirty_price += yield.get_present_value(
                    yearly_coupon / frequency,
                    i / (double)frequency - interval_end,
                    RateType::CONTINUOUS);
        price.dirty_price += yield.get_present_value(principal, n_years - interval_end, RateType::CONTINUOUS);

        const auto lhs_rate = Rate::make_yearly_rate(yearly_coupon / principal, 1);
        price.dirty_price = yield.get_present_value(price.dirty_price, interval_end - point, RateType::CONTINUOUS);
        price.clean_price = price.dirty_price - lhs_rate.get_interest(principal, point - interval_start, RateType::YEARLY);
        return price;
    }

    [[nodiscard]] auto get_yield(int date, double clean_price) const noexcept
    {
        const auto payment_sum = principal + yearly_coupon * n_years;
        if (payment_sum < clean_price)
        {
            std::cerr << "Price of bond can't be more than sum of payments. This implies negative yield which is not reasonable" << std::endl;
            std::terminate();
        }

        double l = 0;
        double r = 1;
        while (std::abs(l - r) > 1e-6)
        {
            const double m = std::midpoint(l, r);
            const double fetched_price = get_price(date, Rate::make_yearly_rate(m , 1)).clean_price;
            if (fetched_price >= clean_price)
                l = m;
            else
                r = m;
        }

        return Rate::make_yearly_rate(l, 1);
    }
};

int main()
{
    std::cout << std::setprecision(8) << std::fixed;

    const auto bond = Bond{
        .principal = 100,
        .yearly_coupon = 3.6,
        .n_years = 15,
        .start_date = 20190430,
        .frequency = 2,
        .convention = DayCountConvention::_30_360
    };

    const auto date = 20241206;
    const auto price = 107;

    std::cout << "Price: " << price << std::endl;
    std::cout << "Yield: \t\t\t" << bond.get_yield(date, price).get_yearly_rate(bond.frequency) * 100 << std::endl;
    std::cout << "Current Yield: \t" << (bond.yearly_coupon / price) * 100 << std::endl;

    return 0;
}
