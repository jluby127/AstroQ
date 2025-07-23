#!/usr/bin/env python3
"""
Example showing how to use the new ScheduleManager instead of data_admin
for the astroq schedule command.
"""

import astroq.schedule_manager as sm
import astroq.splan as splan
import astroq.request as rq


def schedule_with_new_manager(config_path, request_file):
    """
    Example of using the new ScheduleManager for scheduling.
    
    This replaces the old approach:
    manager = mn.data_admin(cf)
    manager.run_admin()
    schedule = sch.Scheduler(request_set, cf)
    """
    
    # Initialize the specialized schedule manager
    schedule_manager = sm.ScheduleManager(config_path)
    
    # Load the request set
    request_set = rq.read_json(request_file)
    
    # Create semester planner
    semester_planner = splan.SemesterPlanner(request_set, config_path)
    
    # Run the model
    semester_planner.run_model()
    
    return semester_planner


def compare_old_vs_new():
    """
    Comparison of old vs new approach.
    """
    
    config_path = "config_2025A.ini"
    request_file = "outputs/2025-07-06/request_set.json"
    
    print("=== OLD APPROACH ===")
    print("""
    # Old way - data_admin does everything
    manager = mn.data_admin(cf)
    manager.run_admin()  # 100+ lines of setup
    schedule = sch.Scheduler(request_set, cf)
    schedule.run_model()
    """)
    
    print("=== NEW APPROACH ===")
    print("""
    # New way - specialized manager
    schedule_manager = sm.ScheduleManager(config_path)
    request_set = rq.read_json(request_file)
    scheduler = sch.Scheduler(request_set, config_path, schedule_manager)
    scheduler.run_model()
    """)
    
    print("=== BENEFITS ===")
    print("1. Focused responsibility - only handles scheduling")
    print("2. Better error handling - validates files exist")
    print("3. Environment variable support - no hardcoded paths")
    print("4. Easier testing - can mock individual components")
    print("5. Clearer interface - explicit methods for getting data")


if __name__ == "__main__":
    compare_old_vs_new() 