"""The Mastermind engine class.

It implements the game logic for the Mastermind game. It handles the generation of
the secret sequence, evaluation of player guesses, and game state management.
"""

import numpy as np
import numpy.typing as npt


class Mastermind:
    """A class to represent the Mastermind game logic."""

    def __init__(self, n: int, k: int, seed: int | None = None):
        """Initialize the game with a sequence length, range of numbers, and optional seed.

        Args:
            n (int): The length of the sequence.
            k (int): The range of numbers (0 to k-1).
            seed (int | None): The seed for random number generation (optional).

        """
        self.seed = seed
        if self.seed is not None:
            np.random.seed(self.seed)
        self.n = n
        self.k = k
        self.restart()
        return

    def restart(self):
        """Generate a new random sequence of numbers based on current n and k."""
        self.tries = 0
        self.secret_sequence = np.random.randint(0, self.k, size=self.n, dtype=np.uint8)
        return

    def evaluate_guess(self, guess: npt.NDArray[np.uint8]) -> tuple[int, int]:
        """Evaluate the player's guess.

        Args:
            guess (npt.NDArray[np.uint8]): The player's guess as a numpy array of uint8.

        Returns:
            int: The number of correct numbers in the correct position.
            int: The number of correct numbers in the wrong position.

        """
        self.tries += 1
        correct_position = (guess == self.secret_sequence).sum()
        correct_number = np.minimum(
            np.bincount(guess, minlength=self.k),
            np.bincount(self.secret_sequence, minlength=self.k),
        ).sum()
        correct_number -= correct_position
        return correct_position, correct_number
