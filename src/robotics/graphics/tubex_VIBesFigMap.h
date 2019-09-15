/** 
 *  \file
 *  VIBesFigMap class
 * ----------------------------------------------------------------------------
 *  \date       2015
 *  \author     Simon Rohou
 *  \copyright  Copyright 2019 Simon Rohou
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#ifndef __TUBEX_VIBESFIGMAP_H__
#define __TUBEX_VIBESFIGMAP_H__

#include <map>
#include <vector>

#include "tubex_VIBesFig.h"
#include "tubex_TubeVector.h"
#include "tubex_TrajectoryVector.h"
#include "tubex_Beacon.h"
#include "tubex_ColorMap.h"

namespace tubex
{
  #define TRAJMAP_NB_DISPLAYED_POINTS   1000
  #define MAP_SLICES_NUMBER_TO_DISPLAY  300000

  // HTML color codes:
  #define DEFAULT_BEACON_COLOR      "#FF5D00[white]"
  #define DEFAULT_TRAJMAP_COLOR     "#276279"
  #define DEFAULT_MAPBCKGRND_COLOR  "#d2d2d2[#d2d2d2]"
  #define DEFAULT_OBS_COLOR         "gray"

  /**
   * \class VIBesFigMap
   * \brief Two-dimensional graphical item to project dynamical items
   *        (tubes, trajectories, etc.) on a map.
   *
   * One figure is linked to some tube or trajectory pointers, so that
   * any update on these objects can be easily displayed on the figure. 
   */
  class VIBesFigMap : public VIBesFig
  {
    public:

      /// \name Basics
      /// @{

      /**
       * \brief Creates a VIBesFigMap
       *
       * \param fig_name name of the figure as displayed in the window title
       */
      VIBesFigMap(const std::string& fig_name);

      /**
       * \brief VIBesFigMap destructor
       */
      ~VIBesFigMap();

      /**
       * \brief Restricts the display of the dynamical items to a part
       *        of their temporal domain only
       *
       * \param restricted_tdomain subset of the temporal domain of the referenced items
       */
      void set_restricted_tdomain(const ibex::Interval& restricted_tdomain);

      /**
       * \brief Enables the display of previous versions of the tubes,
       *        in order to highlight the contractions of these sets
       *
       * \param enable boolean parameter, `true` by default
       */
      void enable_tubes_backgrounds(bool enable = true);

      /**
       * \brief Displays this figure
       */
      void show();

      /**
       * \brief Displays this figure with a custom robot size
       *
       * \param robot_size length of the displayed robot
       */
      void show(float robot_size);
      
      /// @}
      /// \name Handling tubes
      /// @{

      /**
       * \brief Adds a tube vector object to this figure
       *
       * \param tube a const pointer to a TubeVector object to be displayed
       * \param name a name to identify this object
       * \param index_x integer reference for the x-horizontal state component
       * \param index_y integer reference for the y-horizontal state component
       */
      void add_tube(const TubeVector *tube, const std::string& name, int index_x, int index_y);

      /**
       * \brief Sets a new name for a tube
       *
       * \param tube the const pointer to the TubeVector object to be renamed
       * \param name a new name to identify this object
       */
      void set_tube_name(const TubeVector *tube, const std::string& name);

      /**
       * \brief Sets constant color properties for a given tube
       *
       * \param tube the const pointer to the TubeVector object for which the color will be set
       * \param color a color to draw this tube
       */
      void set_tube_color(const TubeVector *tube, const std::string& color);

      /**
       * \brief Sets a color map for a given tube
       *
       * \param tube the const pointer to the TubeVector object for which the color will be set
       * \param colormap a ColorMap to draw this tube with time evolving colors
       * \param traj_colormap an optional const pointer to a Trajectory object in order to have a custom
       *        color evolution instead of a linear time evolving color 
       */
      void set_tube_color(const TubeVector *tube, const ColorMap& colormap, const Trajectory *traj_colormap = NULL);

      /**
       * \brief Removes a tube from this figure
       *
       * \todo automatically clear the figure (for now, the tube is only removed
       *       from the list of objects to be shown afterwards)
       *
       * \note the object will not be deleted
       *
       * \param tube the const pointer to the TubeVector object to be removed
       */
      void remove_tube(const TubeVector *tube);

      /// @}
      /// \name Handling trajectories
      /// @{

      /**
       * \brief Adds a trajectory vector object (x,y) to this figure
       *
       * The heading information may be required for vehicle display. Here it will be simulated from
       * the two components of the traj object. If available, use the other add_trajectory method.
       *
       * \param traj a const pointer to a TrajectoryVector object to be displayed
       * \param name a name to identify this object
       * \param index_x integer reference for the x-horizontal state component
       * \param index_y integer reference for the y-horizontal state component
       * \param color an optional constant color for this object
       */
      void add_trajectory(const TrajectoryVector *traj, const std::string& name, int index_x, int index_y, const std::string& color = DEFAULT_TRAJMAP_COLOR);

      /**
       * \brief Adds a trajectory vector object (x,y,heading) to this figure
       *
       * The heading information may be required for vehicle display. 
       * If not available, use the other add_trajectory method.
       *
       * \param traj a const pointer to a TrajectoryVector object to be displayed
       * \param name a name to identify this object
       * \param index_x integer reference for the x-horizontal state component
       * \param index_y integer reference for the y-horizontal state component
       * \param index_heading integer reference for the heading state component
       * \param color an optional constant color for this object
       */
      void add_trajectory(const TrajectoryVector *traj, const std::string& name, int index_x, int index_y, int index_heading, const std::string& color = DEFAULT_TRAJMAP_COLOR);

      /**
       * \brief Sets a new name for a trajectory
       *
       * \param traj the const pointer to the TrajectoryVector object to be renamed
       * \param name a new name to identify this object
       */
      void set_trajectory_name(const TrajectoryVector *traj, const std::string& name);

      /**
       * \brief Sets constant color properties for a given trajectory
       *
       * \param traj the const pointer to the TrajectoryVector object for which the color will be set
       * \param color a color to draw this trajectory
       */
      void set_trajectory_color(const TrajectoryVector *traj, const std::string& color);

      /**
       * \brief Sets a color map for a given trajectory
       *
       * \param traj the const pointer to the TrajectoryVector object for which the color will be set
       * \param colormap a ColorMap to draw this trajectory with time evolving colors
       * \param traj_colormap an optional const pointer to a Trajectory object in order to have a custom
       *        color evolution instead of a linear time evolving color 
       */
      void set_trajectory_color(const TrajectoryVector *traj, const ColorMap& colormap, const Trajectory *traj_colormap = NULL);

      /**
       * \brief Removes a trajectory from this figure
       *
       * \todo automatically clear the figure (for now, the trajectory is only removed
       *       from the list of objects to be shown afterwards)
       *
       * \note the object will not be deleted
       *
       * \param traj the const pointer to the TrajectoryVector object to be removed
       */
      void remove_trajectory(const TrajectoryVector *traj);

      /// @}
      /// \name Handling robotics objects
      /// @{

      /**
       * \brief Draws a vehicle on top of its trajectory
       *
       * \param t temporal key of the vehicle state to be drawn
       * \param traj the const pointer to the related TrajectoryVector object
       * \param params VIBes parameters related to the vehicle (for groups)
       * \param size optional robot size (-1 = size of main vehicle by default)
       */
      void draw_vehicle(double t, const TrajectoryVector *traj, float size = -1);

      /**
       * \brief Draws a vehicle on top of its trajectory
       *
       * \param t temporal key of the vehicle state to be drawn
       * \param traj the const pointer to the related TrajectoryVector object
       * \param params VIBes parameters related to the vehicle (for groups)
       * \param size optional robot size (-1 = size of main vehicle by default)
       */
      void draw_vehicle(double t, const TrajectoryVector *traj, const vibes::Params& params, float size = -1);

      /**
       * \brief Adds a Beacon object to the map
       *
       * The width of the displayed box is related to the uncertainty of the Beacon object.
       *
       * \param beacon a const reference to the Beacon object to be drawn
       * \param color optional color of the beacon
       */
      void add_beacon(const Beacon& beacon, const std::string& color = DEFAULT_BEACON_COLOR);

      /**
       * \brief Adds a Beacon object to the map with a specific width
       *
       * \param beacon a const reference to the Beacon object to be drawn
       * \param width the real width of the squared beacon
       * \param color optional color of the beacon
       */
      void add_beacon(const Beacon& beacon, double width, const std::string& color = DEFAULT_BEACON_COLOR);

      /**
       * \brief Adds a set of Beacon objects to the map
       *
       * The width of the displayed boxes is related to the uncertainty of the Beacon objects.
       *
       * \param v_beacons a vector of const references to Beacon objects to be drawn
       * \param color optional color of the beacons
       */
      void add_beacons(const std::vector<Beacon>& v_beacons, const std::string& color = DEFAULT_BEACON_COLOR);

      /**
       * \brief Adds a set of Beacon objects to the map with a specific width
       *
       * \param v_beacons a vector of const references to Beacon objects to be drawn
       * \param width the real width of the squared beacons
       * \param color optional color of the beacons
       */
      void add_beacons(const std::vector<Beacon>& v_beacons, double width, const std::string& color = DEFAULT_BEACON_COLOR);
      
      /**
       * \brief Adds a range-and-bearing uncertain observation to the map
       *
       * The observation is a 3d interval vector: (time,range,bearing).
       * It is associated to a trajectory, from which the observation has been made.
       *
       * \param obs the 3d interval vector enclosing the measurement
       * \param traj the const pointer to the related TrajectoryVector object
       * \param color optional color of the observation
       */
      void add_observation(const ibex::IntervalVector& obs, const TrajectoryVector *traj, const std::string& color = DEFAULT_OBS_COLOR);
      
      /**
       * \brief Adds a set of range-and-bearing uncertain observations to the map
       *
       * The observation is a 3d interval vector: (time,range,bearing).
       * It is associated to a trajectory, from which the observation has been made.
       *
       * \param v_obs a vector of 3d interval vectors enclosing the measurements
       * \param traj the const pointer to the related TrajectoryVector object
       * \param color optional color of the observation
       */
      void add_observations(const std::vector<ibex::IntervalVector>& v_obs, const TrajectoryVector *traj, const std::string& color = DEFAULT_OBS_COLOR);
  
      /// @}

    protected:

      /**
       * \brief Draws a tube
       *
       * \param tube the const pointer to the TubeVector object to be drawn
       * \return the projected boxed envelope of the tube (used for figure calibration)
       */
      const ibex::IntervalVector draw_tube(const TubeVector *tube);

      /**
       * \brief Draws the background of a tube
       *
       * Backgrounds are previous versions of the tubes, before their last contraction.
       * This is useful to highlight on the map the improvement provided by the contractions.
       *
       * \param tube the const pointer to the TubeVector object to be drawn
       * \return the projected boxed envelope of the tube background (used for figure calibration)
       */
      const ibex::IntervalVector draw_tube_background(const TubeVector *tube);

      /**
       * \brief Draws a trajectory
       *
       * \todo ColorMap associated to point display (for now: constant color for points)
       *
       * \param traj the const pointer to the TrajectoryVector object to be drawn
       * \param points_size optional size of points, 0 by default (line display instead of a list of points)
       * \return the projected boxed envelope of the trajectory (used for figure calibration)
       */
      const ibex::IntervalVector draw_trajectory(const TrajectoryVector *traj, float points_size = 0.);

      /**
       * \brief Draws the slices of a tube
       *
       * \param tube the const pointer to the TubeVector object to be drawn
       */
      void draw_slices(const TubeVector *tube);

      /**
       * \brief Draws a Beacon object
       *
       * The width of the displayed box is related to the uncertainty of the Beacon object.
       *
       * \param beacon a const reference to the Beacon object to be drawn
       * \param color color of the beacon
       * \param params VIBes parameters related to the beacon
       */
      void draw_beacon(const Beacon& beacon, const std::string& color, const vibes::Params& params);

      /**
       * \brief Draws a Beacon object with a specific width
       *
       * \param beacon a const reference to the Beacon object to be drawn
       * \param width the real width of the squared beacon
       * \param color color of the beacon
       * \param params VIBes parameters related to the beacon
       */
      void draw_beacon(const Beacon& beacon, double width, const std::string& color, const vibes::Params& params);

      /**
       * \brief Draws a range-and-bearing uncertain observation on the map
       *
       * The observation is a 3d interval vector: (time,range,bearing).
       * It is associated to a trajectory, from which the observation has been made.
       *
       * \param obs the 3d interval vector enclosing the measurement
       * \param traj the const pointer to the related TrajectoryVector object
       * \param color color of the observation
       * \param params VIBes parameters related to the observation
       */
      void draw_observation(const ibex::IntervalVector& obs, const TrajectoryVector *traj, const std::string& color, const vibes::Params& params);

    protected:

      /**
       * \struct FigMapTubeParams
       * \brief Specifies some parameters related to a Tube display
       */
      struct FigMapTubeParams
      {
        std::string name; //!< human readable identifier of the tube
        std::string color = ""; //!< constant color of the tube
        std::pair<ColorMap,const Trajectory*> color_map; //!< custom colormap of the tube
        const Tube *tube_x_copy = NULL, *tube_y_copy; //!< to display previous values in background, before any new contraction
        bool from_first_to_last = !true; //!< order of slices display (each one on top of the next/previous one)
        int index_x, index_y; //!< integer references for state components
      };

      /**
       * \struct FigMapTrajParams
       * \brief Specifies some parameters related to a Trajectory display
       */
      struct FigMapTrajParams
      {
        std::string name; //!< human readable identifier of the tube
        std::string color = ""; //!< constant color of the trajectory
        std::pair<ColorMap,const Trajectory*> color_map; //!< custom colormap of the trajectory
        int index_x, index_y, index_heading = -1; //!< integer references for state components ; -1: not defined
      };

      std::map<const TubeVector*,FigMapTubeParams> m_map_tubes; //!< map of Tube objects to be displayed, together with parameters
      std::map<const TrajectoryVector*,FigMapTrajParams> m_map_trajs; //!< map of Trajectory objects to be displayed, together with parameters

      ibex::Interval m_restricted_tdomain; //!< restricts the display to a part of the temporal domain
      bool m_draw_tubes_backgrounds = true; //!< if `true`, will highlight tubes contractions
      float m_robot_size = 5.5; //!< if `0`, no robot display
  };
}

#endif