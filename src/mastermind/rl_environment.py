"""A reinforcement learning environment wrapper for the Mastermind game.

It defines the state space, action space, and reward structure for RL agents.
"""

import numpy as np

from .engine import Mastermind


class MastermindEnv:
    """RL environment wrapper for Mastermind game."""

    def __init__(self, n: int, k: int, seed: int | None = None) -> None:
        """Initialize the RL environment.

        Args:
            n (int): Length of the sequence
            k (int): Range of numbers (0 to k-1)
            seed (int | None): Random seed for reproducibility

        """
        self.game = Mastermind(n, k, seed)
        self.current_guess = []
        self.history = []  # List of (guess, feedback) tuples
        return

    def step(self, action: int) -> tuple[dict, float, bool]:
        """Take a step in the environment by adding a number to the current guess.

        Args:
            action (int): The number to add (0 to k-1)

        Returns:
            dict: Current state observation
            float: Reward (-1, 0, or +1)
            bool: Whether the episode is done

        """
        self.current_guess.append(action)

        # If guess is not complete yet
        if len(self.current_guess) < self.game.n:
            return self._get_state(), 0.0, False

        # Complete guess, evaluate it
        guess_array = np.array(self.current_guess, dtype=np.uint8)
        correct_pos, correct_num = self.game.evaluate_guess(guess_array)

        # Store in history
        self.history.append((guess_array, (correct_pos, correct_num)))

        # Determine reward
        reward = 1.0 if correct_pos == self.game.n else -1.0
        done = correct_pos == self.game.n

        # Reset current guess
        self.current_guess = []

        return self._get_state(), reward, done

    def reset(self) -> dict:
        """Reset the environment for a new episode.

        Returns:
            dict: Initial state observation

        """
        self.game.restart()
        self.current_guess = []
        self.history = []
        return self._get_state()

    def _get_state(self) -> dict:
        """Get the current state observation.

        Returns:
            dict: State containing current position, history, and available actions

        """
        return {
            "position": len(self.current_guess),
            "history": self.history,
        }
