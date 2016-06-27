//
// Created by clonker on 07.03.16.
//

#include <boost/log/core/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/thread/thread.hpp>
#include <readdy/common/Utils.h>
#include "gtest/gtest.h"

int perform_tests(int argc, char **argv) {
    readdy::utils::testing::loadPlugins();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

int main(int argc, char **argv) {
    namespace expr = boost::log::expressions;
    namespace attr = boost::log::attributes;
    namespace dam = boost::log::aux::default_attribute_names;
    namespace keywords = boost::log::keywords;

    auto coreLogger = boost::log::core::get();
    coreLogger->add_global_attribute(dam::timestamp(), attr::local_clock());
    coreLogger->add_global_attribute(dam::thread_id(), attr::current_thread_id());
    auto loggingStream = expr::stream
                         << "[          ] "
                         << "[" << expr::format_date_time<boost::posix_time::ptime>("TimeStamp", "%Y-%m-%d %H:%M:%S") << "] "
                         << "[" << expr::attr<attr::current_thread_id::value_type>("ThreadID") << "] "
                         << "[" << boost::log::trivial::severity << "] "
                         << expr::smessage;
    boost::log::add_console_log(std::cout, keywords::format = loggingStream);

    int result = perform_tests(argc, argv);
    boost::log::core::get()->flush();
    boost::log::core::get()->remove_all_sinks();
    return result;
}
