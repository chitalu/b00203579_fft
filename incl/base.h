#ifndef __BASE_H__
#define __BASE_H__


struct cpu_timer_t {
  typedef std::chrono::time_point<std::chrono::system_clock,
                                  std::chrono::system_clock::duration>
      cpu_prof_time_t;

  cpu_timer_t(std::string prof_instance) : m_prof_instance(prof_instance) {
    m_start_time = std::chrono::high_resolution_clock::now();
  }

  ~cpu_timer_t(void) {
    m_end_time = std::chrono::high_resolution_clock::now();

    profile_time_t elapsed =
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            m_end_time - m_start_time).count();
    g_profiling_data[m_prof_instance].push_back(elapsed);

    if (g_profiling_data[m_prof_instance].size() > g_profiling_samples) {
      g_profiling_data[m_prof_instance].pop_front();
    }
  }

private:
  std::string m_prof_instance;
  cpu_prof_time_t m_start_time, m_end_time;
};

#define PROFILE_()                                                    		   \
  std::shared_ptr<cpu_timer_t> profiler_;                                 \
  if (g_profiling_enabled) {                                                   \
    std::shared_ptr<cpu_timer_t> p_inst =                                 \
        std::shared_ptr<cpu_timer_t>(new cpu_timer_t(__FUNCTION__)); \
    profiler_ = p_inst;                                                        \
  }


#endif