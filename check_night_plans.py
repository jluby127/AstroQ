#!/usr/bin/env python3
"""
Check night_plan.csv files in band directories.

This script opens a series of folders (one for each band) at the holders_dir
location and checks:
1. If night_plan.csv exists
2. Counts how many stars it has
"""

import os
import sys
import argparse
from datetime import datetime
import pandas as pd
from pathlib import Path


def extract_log_statistics(log_path):
    """Extract statistics from the log file.
    
    Args:
        log_path: Path to the log file
        
    Returns:
        dict: Dictionary with extracted statistics or None if not found
    """
    if not log_path.exists():
        return None
    
    stats = {}
    
    try:
        with open(log_path, 'r') as f:
            content = f.read()
            
            # Extract requested vs scheduled
            import re
            requested_match = re.search(r'Observations Requested:\s+(\d+)', content)
            scheduled_match = re.search(r'Observations Scheduled:\s+(\d+)', content)
            
            if requested_match:
                stats['requested'] = int(requested_match.group(1))
            if scheduled_match:
                stats['scheduled'] = int(scheduled_match.group(1))
            
            # Extract observing duration
            duration_match = re.search(r'Observing Duration \(min\):\s+([\d.]+)', content)
            if duration_match:
                stats['duration_min'] = float(duration_match.group(1))
            
            # Extract time spent
            exposing_match = re.search(r'Time Spent Exposing \(min\):\s+([\d.]+)', content)
            if exposing_match:
                stats['exposing_min'] = float(exposing_match.group(1))
            
            idle_match = re.search(r'Time Spent Idle \(min\):\s+([\d.]+)', content)
            if idle_match:
                stats['idle_min'] = float(idle_match.group(1))
            
            slew_match = re.search(r'Time Spent Slewing \(min\):\s+([\d.]+)', content)
            if slew_match:
                stats['slew_min'] = float(slew_match.group(1))
            
            # Extract total open shutter time
            shutter_match = re.search(r'Total Open Shutter Time Scheduled:\s+([\d.]+)\s+hours', content)
            if shutter_match:
                stats['shutter_hours'] = float(shutter_match.group(1))
                
    except Exception as e:
        print(f"Warning: Error reading log file {log_path}: {e}")
    
    return stats if stats else None


def check_night_plan(base_path, band, output_path=None, date=None):
    """Check if night_plan.csv exists and count the number of stars.
    
    Args:
        base_path: Base directory path for CSV files (should include semester/date)
        band: Band name (e.g., 'band1', 'band2', etc.)
        output_path: Base directory path for log files (CC_OUTPUT_PATH/semester/date)
        date: Date string in YYYY-MM-DD format for checking allocation.csv
    
    Returns:
        tuple: (exists: bool, star_count: int, log_stats: dict, error_message: str)
    """
    band_path = Path(base_path) / band
    csv_path = band_path / "output" / "night_plan.csv"
    
    # Log path is in the output directory structure
    if output_path:
        log_path = Path(output_path) / band / "astroq.log"
    else:
        log_path = band_path / "astroq.log"
    
    # Check if band directory exists
    if not band_path.exists():
        return False, 0, None, f"Band directory does not exist"
    
    # Check if output directory exists
    output_dir = band_path / "output"
    if not output_dir.exists():
        return False, 0, None, f"Output directory does not exist"
    
    # Check if night_plan.csv exists
    if not csv_path.exists():
        # Before reporting error, check if date is in allocation.csv
        if date:
            # Look for allocation.csv in the band directory or parent directory
            allocation_paths = [
                band_path / "allocation.csv",
                base_path / "allocation.csv",
                band_path.parent / "allocation.csv",
            ]
            
            date_in_allocation = False
            for alloc_path in allocation_paths:
                if alloc_path.exists():
                    try:
                        alloc_df = pd.read_csv(alloc_path)
                        # Check if 'start' column exists (as used in nplan.py)
                        if 'start' in alloc_df.columns:
                            # Extract date part from start column (format: YYYY-MM-DD or YYYY-MM-DD HH:MM:SS)
                            date_str = str(date)
                            # Check if date matches any start date (compare first 10 chars for YYYY-MM-DD)
                            alloc_df['start_date'] = alloc_df['start'].astype(str).str[:10]
                            if date_str in alloc_df['start_date'].values:
                                date_in_allocation = True
                                break
                        # Also check for other common date column names
                        for col in alloc_df.columns:
                            if col.lower() in ['date', 'night']:
                                date_str = str(date)
                                if date_str in alloc_df[col].astype(str).values:
                                    date_in_allocation = True
                                    break
                        if date_in_allocation:
                            break
                    except Exception as e:
                        # If we can't read the file, continue checking other paths
                        continue
            
            # If date is not in allocation.csv, this is expected - no CC time tonight
            if not date_in_allocation:
                return False, 0, None, "No CC Time Tonight"
        
        # If we can't check allocation or date is in allocation but file still missing, report error
        return False, 0, None, f"night_plan.csv does not exist"
    
    # Try to read and count stars
    star_count = 0
    try:
        df = pd.read_csv(csv_path)
        star_count = len(df)
    except Exception as e:
        return False, 0, None, f"Error reading CSV: {str(e)}"
    
    # Extract log statistics
    log_stats = extract_log_statistics(log_path)
    
    return True, star_count, log_stats, ""


def main():
    """Main function to check all band directories."""
    
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Check night_plan.csv files in band directories',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '-s', '--semester',
        type=str,
        default='2025B',
        help='Semester string (e.g., 2025B). Default: 2025B'
    )
    parser.add_argument(
        '-d', '--date',
        type=str,
        default=None,
        help='Date string in YYYY-MM-DD format (e.g., 2025-10-29). Default: today\'s date'
    )
    parser.add_argument(
        'holders_dir',
        nargs='?',
        default=None,
        help='Optional: Full path to holders directory. If not provided, will construct from environment variables and flags.'
    )
    
    args = parser.parse_args()
    
    # Default bands list
    bands = ['band1', 'band2', 'band3', 'full-band1', 'full-band2', 'full-band3']
    
    # Get semester and date from arguments
    semester = args.semester
    
    # Get date - use provided date or default to today
    if args.date is None:
        date = datetime.now().strftime('%Y-%m-%d')
    else:
        date = args.date
        # Validate date format
        try:
            datetime.strptime(date, '%Y-%m-%d')
        except ValueError:
            print(f"❌ Error: Invalid date format '{date}'. Please use YYYY-MM-DD format (e.g., 2025-10-29)")
            sys.exit(1)
    
    # Get holders directory
    if args.holders_dir:
        holders_dir = Path(args.holders_dir)
    else:
        # Try to construct from environment variables
        if 'CC_RESULTS_COPY_PATH' in os.environ:
            base_path = os.environ['CC_RESULTS_COPY_PATH']
            holders_dir = Path(base_path) / semester / date
        else:
            print("Usage: python check_night_plans.py [-s SEMESTER] [-d DATE] [holders_dir]")
            print("\nOptions:")
            print("  -s, --semester   Semester string (e.g., 2025B). Default: 2025B")
            print("  -d, --date       Date string in YYYY-MM-DD format. Default: today's date")
            print("  holders_dir      Optional: Full path to holders directory")
            print("\nOr set environment variables:")
            print("  CC_RESULTS_COPY_PATH - Base path for results")
            sys.exit(1)
    
    holders_dir = Path(holders_dir)
    
    if not holders_dir.exists():
        print(f"❌ Error: Holders directory does not exist: {holders_dir}")
        print(f"\nConstructed path from:")
        if 'CC_RESULTS_COPY_PATH' in os.environ:
            print(f"  CC_RESULTS_COPY_PATH: {os.environ['CC_RESULTS_COPY_PATH']}")
        print(f"  SEMESTER: {semester}")
        print(f"  DATE: {date}")
        sys.exit(1)
    
    # Construct output path for log files using the same semester and date
    output_dir = None
    if 'CC_OUTPUT_PATH' in os.environ:
        output_base = Path(os.environ['CC_OUTPUT_PATH'])
        output_dir = output_base / semester / date
        if output_dir.exists():
            print(f"📊 Checking logs in: {output_dir}")
    
    print(f"📁 Checking night_plan.csv files in: {holders_dir}\n")
    
    results = []
    
    # Check each band
    for band in bands:
        exists, star_count, log_stats, error = check_night_plan(holders_dir, band, output_dir, date)
        results.append((band, exists, star_count, log_stats, error))
    
    # Write results to file
    output_file = holders_dir / "quicklook_status.txt"
    
    with open(output_file, 'w') as f:
        f.write("=" * 145 + "\n")
        # Multi-line header for better alignment
        f.write(f"{'':<15} {'':<8} {'N':<12} {'N':<12} {'Time in':<15} {'Time':<12} {'Time':<12} {'Time':<12} {'Open':<15} {'':<20}\n")
        f.write(f"{'':<15} {'':<8} {'Stars':<12} {'Stars':<12} {'Night':<15} {'Exposing':<12} {'Slewing':<12} {'Idle':<12} {'Shutter':<15} {'':<20}\n")
        f.write(f"{'Band':<15} {'Exists':<8} {'Requested':<12} {'Scheduled':<12} {'(min)':<15} {'(min)':<12} {'(min)':<12} {'(min)':<12} {'(hrs)':<15} {'Status':<20}\n")
        f.write("=" * 145 + "\n")
        
        # First pass: determine status for each band and collect for summary
        band_statuses = []
        for band, exists, star_count, log_stats, error in results:
            if exists:
                status = "OK"
            elif error == "No CC Time Tonight":
                status = "N/A"
            else:
                status = f"ERROR: {error}"
            
            # Format log statistics
            if log_stats:
                requested = log_stats.get('requested', 'N/A')
                scheduled = log_stats.get('scheduled', 'N/A')
                duration = f"{log_stats.get('duration_min', 'N/A'):.1f}" if 'duration_min' in log_stats else "N/A"
                exposing = f"{log_stats.get('exposing_min', 'N/A'):.1f}" if 'exposing_min' in log_stats else "N/A"
                slew = f"{log_stats.get('slew_min', 'N/A'):.1f}" if 'slew_min' in log_stats else "N/A"
                idle = f"{log_stats.get('idle_min', 'N/A'):.1f}" if 'idle_min' in log_stats else "N/A"
                shutter = f"{log_stats.get('shutter_hours', 'N/A'):.2f}" if 'shutter_hours' in log_stats else "N/A"
                
                # Check if idle time exceeds threshold # of minutes and add warning to status
                warn_idle = 60
                if exists and 'idle_min' in log_stats:
                    idle_min = log_stats.get('idle_min')
                    if idle_min > warn_idle:
                        status = f"WARNING: Idle > {warn_idle}min"
            else:
                requested = "N/A"
                scheduled = "N/A"
                duration = "N/A"
                exposing = "N/A"
                slew = "N/A"
                idle = "N/A"
                shutter = "N/A"
            
            # Store status for summary calculation
            band_statuses.append((exists, status))
            f.write(f"{band:<15} {str(exists):<8} {str(requested):<12} {str(scheduled):<12} {duration:<15} {exposing:<12} {slew:<12} {idle:<12} {shutter:<15} {status:<20}\n")
        
        f.write("=" * 200 + "\n")
        
        # Summary - count based on final status strings
        total_bands = len(results)
        # Successful: bands that exist (regardless of warnings)
        successful = sum(1 for exists, status in band_statuses if exists)
        # Failures: only actual errors (status starting with "ERROR:")
        failures = sum(1 for exists, status in band_statuses if status.startswith("ERROR:"))
        # N/A: bands with "N/A" status (no CC time tonight)
        na_count = sum(1 for exists, status in band_statuses if status == "N/A")
        # Warnings: warnings that aren't "N/A" (e.g., "WARNING: Idle > 60min")
        warnings = sum(1 for exists, status in band_statuses if status.startswith("WARNING:"))
        
        f.write(f"\nSummary:\n")
        f.write(f"   Successful: {successful}\n")
        f.write(f"   Failures: {failures}\n")
        f.write(f"   N/A: {na_count}\n")
        f.write(f"-------------------\n")
        f.write(f"   Total bands: {total_bands}\n")
        f.write(f"   Warnings: {warnings}\n")

    
    print(f"✅ Results written to: {output_file}")
    
    # Always exit with success code - errors are reported in quicklook_status.txt
    # This prevents the script from failing the make command while still documenting issues
    sys.exit(0)


if __name__ == "__main__":
    main()

