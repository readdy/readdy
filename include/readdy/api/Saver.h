//
// Created by mho on 6/12/19.
//

#pragma once

//#ifdef __cpp_lib_filesystem
//#define CPP_FS 1
//#include <filesystem>
//namespace fs = std::filesystem;
//#else
#define CPP_FS 0
#include <readdy/common/filesystem.h>
namespace fs = readdy::util::fs;
//#endif

#include <memory>
#include <type_traits>
#include <utility>

#include <h5rd/h5rd.h>

#include <readdy/common/common.h>
#include <readdy/model/Kernel.h>
#include <readdy/model/IOUtils.h>

namespace readdy::api {

class Saver {
public:
    Saver(std::string base, std::size_t maxNSaves, std::string checkpointTemplate = "checkpoint_{}.h5")
          : _basePath(std::move(base)), _maxNSaves(maxNSaves), _checkpointTemplate(std::move(checkpointTemplate)) {
        {
            // if template is invalid this will raise
            auto testFormat = fmt::format("{:{}}", checkpointTemplate, 123);
        }
        if(fs::exists(_basePath)) {
            // basePath exists, make sure it is a directory
            if(!fs::is_directory(_basePath)) {
                throw std::invalid_argument(fmt::format("Base path \"{}\" exists but is no directory.", _basePath));
            }
            // else we are o.k.!
        } else {
#if CPP_FS
            // basePath did not exist, create it
            if(!std::filesystem::create_directories(basePath)) {
                throw std::invalid_argument(fmt::format("Could not create directory at \"{}\"", basePath));
            }
#else
            throw std::invalid_argument(fmt::format("Base path \"{}\" did not exist", _basePath));
#endif
        }
    }

    void makeCheckpoint(model::Kernel *const kernel, TimeStep t) {
        auto fileName = fmt::format("{:{}}", _checkpointTemplate, t);
        auto filePath = _basePath + "/" + fileName;

        if (_maxNSaves > 0) {
            previousCheckpoints.push(filePath);
        }
        auto file = File::create(filePath, File::Flag::OVERWRITE);
        {

            // write config into checkpoint
            auto cfgGroup = file->createGroup("readdy/config");
            model::ioutils::writeSimulationSetup(cfgGroup, kernel->context());
        }
        {
            model::observables::FlatTrajectory traj(kernel, 1, false);
            traj.enableWriteToFile(*file, "trajectory_ckpt", 1);
            traj.call(t);
        }
        {
            model::observables::Topologies tops(kernel, 1, false);
            tops.enableWriteToFile(*file, "topologies_ckpt", 1);
            tops.call(t);
        }

        while (_maxNSaves > 0 && previousCheckpoints.size() > _maxNSaves) {
            const auto &oldestCheckpoint = previousCheckpoints.front();
            if (fs::exists(oldestCheckpoint)) {
                if (!fs::remove(oldestCheckpoint)) {
                    throw std::runtime_error(fmt::format("Could not remove checkpoint {}", oldestCheckpoint));
                }
            } else {
                log::warn("Tried removing checkpoint {} but it didn't exist (anymore).", oldestCheckpoint);
            }

            previousCheckpoints.pop();
        }
    }

    [[nodiscard]] std::string basePath() const {
        return _basePath;
    }

    [[nodiscard]] std::size_t maxNSaves() const {
        return _maxNSaves;
    }

    [[nodiscard]] std::string checkpointTemplate() const {
        return _checkpointTemplate;
    }

    std::string describe() const {
        std::string description;
        description += fmt::format("   * base path: {}\n", basePath());
        description += fmt::format("   * checkpoint filename template: {}\n", checkpointTemplate());
        description += fmt::format("   * maximal number saves: {}\n", maxNSaves());
        return description;
    }

private:
    std::string _basePath;
    std::size_t _maxNSaves;
    std::string _checkpointTemplate;
    std::queue<std::string> previousCheckpoints {};
};

}
