//
// Created by mho on 6/12/19.
//

#pragma once

#ifdef __cpp_lib_filesystem
#define CPP_FS 1
#include <filesystem>
namespace fs = std::filesystem;
#else
#define CPP_FS 0
#include <readdy/common/filesystem.h>
namespace fs = readdy::util::fs;
#endif

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
          : basePath(std::move(base)), maxNSaves(maxNSaves), checkpointTemplate(std::move(checkpointTemplate)) {
        {
            // if template is invalid this will raise
            auto testFormat = fmt::format(checkpointTemplate, 123);
        }
        if(fs::exists(basePath)) {
            // basePath exists, make sure it is a directory
            if(!fs::is_directory(basePath)) {
                throw std::invalid_argument(fmt::format("Base path \"{}\" exists but is no directory.", basePath));
            }
            // else we are o.k.!
        } else {
#if CPP_FS
            // basePath did not exist, create it
            if(!std::filesystem::create_directories(basePath)) {
                throw std::invalid_argument(fmt::format("Could not create directory at \"{}\"", basePath));
            }
#else
            throw std::invalid_argument(fmt::format("Base path \"{}\" did not exist", basePath));
#endif
        }
    }

    void makeCheckpoint(model::Kernel *const kernel, time_step_type t) {
        auto fileName = fmt::format(checkpointTemplate, t);
        auto filePath = basePath + "/" + fileName;

        if (maxNSaves > 0) {
            previousCheckpoints.push(filePath);
        }

        auto file = File::create(filePath, File::Flag::OVERWRITE);
        {
            // write config into checkpoint
            auto cfgGroup = file->createGroup("readdy/config");
            model::ioutils::writeSimulationSetup(cfgGroup, kernel->context());
        }
        {
            model::observables::FlatTrajectory traj (kernel, 1);
            traj.enableWriteToFile(*file, "trajectory_ckpt", 1);
            traj.callback(t);
        }
        {
            model::observables::Topologies tops (kernel, 1);
            tops.enableWriteToFile(*file, "topologies_ckpt", 1);
            tops.callback(t);
        }

        while(maxNSaves > 0 && previousCheckpoints.size() > maxNSaves) {
            const auto& oldestCheckpoint = previousCheckpoints.front();

            if (fs::exists(oldestCheckpoint)) {
                if(!fs::remove(oldestCheckpoint)) {
                    throw std::runtime_error(fmt::format("Could not remove checkpoint {}", oldestCheckpoint));
                }
            } else {
                log::warn("Tried removing checkpoint {} but it didn't exist (anymore).", oldestCheckpoint);
            }

            previousCheckpoints.pop();
        }
    }

private:
    std::string basePath;
    std::size_t maxNSaves;
    std::string checkpointTemplate;
    std::queue<std::string> previousCheckpoints {};
};

}
