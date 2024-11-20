import numpy as np

from src.mastermind.engine import Mastermind


def test_secret_sequence_with_seed():
    """Test the secret sequence generation with a specific seed."""
    n = 4
    k = 6
    seed = 42
    game = Mastermind(n, k, seed)

    expected_sequence = np.array([4, 1, 3, 5], dtype=np.uint8)

    assert (
        game.secret_sequence == expected_sequence
    ).all(), f"Expected {expected_sequence}, got {game.secret_sequence}"

    guess = np.array([2, 5, 3, 4], dtype=np.uint8)
    correct, misplaced = game.evaluate_guess(guess)

    assert correct == 1, f"Expected 1 correct number, got {correct}"
    assert misplaced == 2, f"Expected 2 misplaced numbers, got {misplaced}"
