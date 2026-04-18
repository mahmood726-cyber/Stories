/**
 * Learner State Management Module
 * Handles all state persistence, badges, streaks, and adaptive learning
 */

const LearnerState = (function() {
  const STORAGE_KEY = 'evidenceReversalLearnerState';

  // Badge definitions with requirements
  const BADGES = {
    first_steps: {
      name: 'First Steps',
      icon: '👣',
      points: 10,
      description: 'Complete first story',
      requirement: (state) => Object.keys(state.completedStories).length > 0
    },
    reflective_practitioner: {
      name: 'Reflective Practitioner',
      icon: '🪞',
      points: 50,
      description: 'Complete 10 reflection journals',
      requirement: (state) => state.reflections.length >= 10
    },
    epistemic_humility: {
      name: 'Epistemic Humility',
      icon: '🎯',
      points: 50,
      description: 'Acknowledge 5 surprise predictions',
      requirement: (state) => state.surprises >= 5
    },
    empathy_navigator: {
      name: 'Empathy Navigator',
      icon: '🧑‍🦳',
      points: 100,
      description: '10 stories as Patient',
      requirement: (state) => countPersonaCompletions(state, 'patient') >= 10
    },
    clinical_translator: {
      name: 'Clinical Translator',
      icon: '👨‍⚕️',
      points: 100,
      description: '10 stories as Clinician',
      requirement: (state) => countPersonaCompletions(state, 'clinician') >= 10
    },
    methods_detective: {
      name: 'Methods Detective',
      icon: '🔬',
      points: 100,
      description: '10 stories as Researcher',
      requirement: (state) => countPersonaCompletions(state, 'researcher') >= 10
    },
    policy_architect: {
      name: 'Policy Architect',
      icon: '🏛️',
      points: 100,
      description: '10 stories as Regulator',
      requirement: (state) => countPersonaCompletions(state, 'regulator') >= 10
    },
    ethics_guardian: {
      name: 'Ethics Guardian',
      icon: '⚖️',
      points: 100,
      description: '10 stories as Ethicist',
      requirement: (state) => countPersonaCompletions(state, 'ethicist') >= 10
    },
    scenario_solver: {
      name: 'Scenario Solver',
      icon: '🎭',
      points: 150,
      description: 'Complete 5 branching scenarios',
      requirement: (state) => Object.keys(state.scenarioCompleted || {}).filter(k => state.scenarioCompleted[k]).length >= 5
    },
    streak_master: {
      name: 'Streak Master',
      icon: '🔥',
      points: 75,
      description: '7-day streak',
      requirement: (state) => state.streak.longest >= 7
    },
    peer_teacher: {
      name: 'Peer Teacher',
      icon: '👥',
      points: 75,
      description: '5 peer teaching explanations',
      requirement: (state) => state.peerTeachings >= 5
    },
    commitment_keeper: {
      name: 'Commitment Keeper',
      icon: '✅',
      points: 50,
      description: 'Set and report on 3 commitments',
      requirement: (state) => state.commitments.filter(c => c.completed).length >= 3
    },
    renaissance_learner: {
      name: 'Renaissance Learner',
      icon: '🌟',
      points: 500,
      description: 'All 5 personas on 10 stories',
      requirement: (state) => {
        const personas = ['patient', 'clinician', 'researcher', 'regulator', 'ethicist'];
        return personas.every(p => countPersonaCompletions(state, p) >= 10);
      }
    },
    evidence_master: {
      name: 'Evidence Master',
      icon: '👑',
      points: 1000,
      description: 'All badges earned',
      requirement: (state) => {
        const otherBadges = Object.keys(BADGES).filter(b => b !== 'evidence_master');
        return otherBadges.every(b => state.badges.includes(b));
      }
    }
  };

  function countPersonaCompletions(state, persona) {
    let count = 0;
    for (const storyId in state.completedStories) {
      if (state.completedStories[storyId]?.includes(persona)) {
        count++;
      }
    }
    return count;
  }

  function getDefaultState() {
    return {
      completedStories: {}, // { "1": ["patient", "clinician"], "2": ["researcher"] }
      currentPersona: 'patient',
      reflections: [], // [{ storyId, persona, text, timestamp, assumptions: [] }]
      predictions: [], // [{ storyId, predicted, actual, confidence, timestamp }]
      commitments: [], // [{ text, dueDate, completed, storyId, timestamp }]
      streak: {
        current: 0,
        longest: 0,
        lastDate: null,
        freezeAvailable: true,
        freezeUsedThisWeek: false
      },
      points: 0,
      badges: [],
      surprises: 0,
      peerTeachings: 0,
      journalEntries: [],
      adaptiveProfile: {
        strengths: [],
        weaknesses: [],
        preTestCompleted: false,
        preTestResults: {}
      },
      scenarioCompleted: {},
      scenarioPoints: {},
      stageProgress: {}, // { storyId_persona: currentStage }
      lastActivity: null
    };
  }

  function loadState() {
    try {
      const stored = localStorage.getItem(STORAGE_KEY);
      if (stored) {
        const parsed = JSON.parse(stored);
        // Merge with defaults to handle new fields
        return { ...getDefaultState(), ...parsed };
      }
    } catch (e) {
      console.error('Failed to load learner state:', e);
    }
    return getDefaultState();
  }

  function saveState(state) {
    try {
      localStorage.setItem(STORAGE_KEY, JSON.stringify(state));
    } catch (e) {
      console.error('Failed to save learner state:', e);
    }
  }

  // Streak management
  function updateStreak(state) {
    const today = new Date().toDateString();
    const lastDate = state.streak.lastDate;

    if (!lastDate) {
      state.streak.current = 1;
      state.streak.lastDate = today;
    } else if (lastDate === today) {
      // Same day, no change
    } else {
      const last = new Date(lastDate);
      const now = new Date(today);
      const diffDays = Math.floor((now - last) / (1000 * 60 * 60 * 24));

      if (diffDays === 1) {
        // Consecutive day
        state.streak.current++;
      } else if (diffDays === 2 && state.streak.freezeAvailable && !state.streak.freezeUsedThisWeek) {
        // Use freeze
        state.streak.freezeUsedThisWeek = true;
        state.streak.freezeAvailable = false;
      } else {
        // Streak broken
        state.streak.current = 1;
      }
      state.streak.lastDate = today;
    }

    // Update longest streak
    if (state.streak.current > state.streak.longest) {
      state.streak.longest = state.streak.current;
    }

    // Reset freeze weekly
    const weekStart = getWeekStart();
    if (!state.streak.weekStart || state.streak.weekStart !== weekStart) {
      state.streak.weekStart = weekStart;
      state.streak.freezeUsedThisWeek = false;
      state.streak.freezeAvailable = true;
    }

    return state;
  }

  function getWeekStart() {
    const now = new Date();
    const dayOfWeek = now.getDay();
    const diff = now.getDate() - dayOfWeek;
    return new Date(now.setDate(diff)).toDateString();
  }

  // Confidence-weighted scoring
  function calculatePoints(isCorrect, confidence) {
    const basePoints = isCorrect ? 30 : 0;

    if (confidence >= 8) {
      // High confidence
      return isCorrect ? Math.round(basePoints * 1.5) : Math.round(basePoints * 0.5);
    } else if (confidence <= 3) {
      // Low confidence - no penalty for uncertainty
      return basePoints;
    } else {
      // Medium confidence
      return basePoints;
    }
  }

  // Badge checking and awarding
  function checkBadges(state) {
    const newBadges = [];

    for (const [badgeId, badge] of Object.entries(BADGES)) {
      if (!state.badges.includes(badgeId) && badge.requirement(state)) {
        state.badges.push(badgeId);
        state.points += badge.points;
        newBadges.push(badgeId);
      }
    }

    return newBadges;
  }

  // Mark story completed for a persona
  function completeStory(state, storyId, persona) {
    if (!state.completedStories[storyId]) {
      state.completedStories[storyId] = [];
    }
    if (!state.completedStories[storyId].includes(persona)) {
      state.completedStories[storyId].push(persona);
    }
    state.lastActivity = new Date().toISOString();
    return updateStreak(state);
  }

  // Add reflection
  function addReflection(state, storyId, persona, text, assumptions = []) {
    state.reflections.push({
      storyId,
      persona,
      text,
      assumptions,
      timestamp: new Date().toISOString()
    });
    return state;
  }

  // Add prediction
  function addPrediction(state, storyId, predicted, actual, confidence) {
    const isSurprised = predicted !== actual && confidence >= 7;
    if (isSurprised) {
      state.surprises++;
    }

    state.predictions.push({
      storyId,
      predicted,
      actual,
      confidence,
      surprised: isSurprised,
      timestamp: new Date().toISOString()
    });

    // Calculate and add points
    const isCorrect = predicted === actual;
    const pointsEarned = calculatePoints(isCorrect, confidence);
    state.points += pointsEarned;

    return { state, pointsEarned, isSurprised };
  }

  // Add commitment
  function addCommitment(state, storyId, text, dueDate) {
    state.commitments.push({
      storyId,
      text,
      dueDate,
      completed: false,
      timestamp: new Date().toISOString()
    });
    return state;
  }

  // Complete commitment
  function completeCommitment(state, index) {
    if (state.commitments[index]) {
      state.commitments[index].completed = true;
      state.commitments[index].completedAt = new Date().toISOString();
    }
    return state;
  }

  // Add peer teaching
  function addPeerTeaching(state, storyId, persona, explanation) {
    state.peerTeachings++;
    state.lastActivity = new Date().toISOString();
    return state;
  }

  // Journal entry management
  function addJournalEntry(state, entry) {
    state.journalEntries.push({
      ...entry,
      id: Date.now(),
      timestamp: new Date().toISOString()
    });
    return state;
  }

  function deleteJournalEntry(state, entryId) {
    state.journalEntries = state.journalEntries.filter(e => e.id !== entryId);
    return state;
  }

  // Stage progress management
  function getStageProgress(state, storyId, persona) {
    const key = `${storyId}_${persona}`;
    return state.stageProgress[key] || 1;
  }

  function setStageProgress(state, storyId, persona, stage) {
    const key = `${storyId}_${persona}`;
    state.stageProgress[key] = stage;
    return state;
  }

  // Export/Import functionality
  function exportState(state) {
    return JSON.stringify(state, null, 2);
  }

  function importState(jsonString) {
    try {
      const imported = JSON.parse(jsonString);
      const merged = { ...getDefaultState(), ...imported };
      saveState(merged);
      return { success: true, state: merged };
    } catch (e) {
      return { success: false, error: e.message };
    }
  }

  // Statistics
  function getStatistics(state) {
    const totalStories = 40;
    const completedCount = Object.keys(state.completedStories).length;
    const personas = ['patient', 'clinician', 'researcher', 'regulator', 'ethicist'];
    const personaCounts = {};
    personas.forEach(p => {
      personaCounts[p] = countPersonaCompletions(state, p);
    });

    const totalBadges = Object.keys(BADGES).length;
    const earnedBadges = state.badges.length;

    const todayDate = new Date().toDateString();
    const todayJournalEntries = state.journalEntries.filter(
      e => new Date(e.timestamp).toDateString() === todayDate
    ).length;

    return {
      completedCount,
      totalStories,
      personaCounts,
      streak: state.streak,
      points: state.points,
      badges: earnedBadges,
      totalBadges,
      journalTotal: state.journalEntries.length,
      journalToday: todayJournalEntries,
      surprises: state.surprises,
      commitments: state.commitments.length,
      completedCommitments: state.commitments.filter(c => c.completed).length
    };
  }

  // Public API
  return {
    BADGES,
    loadState,
    saveState,
    updateStreak,
    calculatePoints,
    checkBadges,
    completeStory,
    addReflection,
    addPrediction,
    addCommitment,
    completeCommitment,
    addPeerTeaching,
    addJournalEntry,
    deleteJournalEntry,
    getStageProgress,
    setStageProgress,
    exportState,
    importState,
    getStatistics,
    countPersonaCompletions
  };
})();

// Export for use in main file
if (typeof window !== 'undefined') {
  window.LearnerState = LearnerState;
}
