#pragma once

#include <vector>
#include <map>
#include <string>
#include <unordered_map>
#include <memory>

// For profiling time spent in a particular function
#define ProfilerPush() Profile::Push(__FUNCTION__)
#define ProfilerPop() Profile::Pop(__FUNCTION__)



class Profiler {
public:
    typedef std::shared_ptr<Profiler> Ptr_t;

    virtual void Push(const std::string& name);
    virtual void Pop(const std::string& name);

    virtual bool PrintTo(const std::string& filename);
    virtual void Print();
};




namespace Profile {
    void Push(const std::string& name);
    void Pop(const std::string& name);

    bool PrintTo(const std::string& filename);
    void Print();

    void SetProfiler(Profiler* profiler);
};