import math

def align(
    seq1: str,
    seq2: str,
    match_award=-3,
    indel_penalty=5,
    sub_penalty=1,
    banded_width=-1,
    gap='-'
) -> tuple[float, str | None, str | None]:

    def compare_chars(char1, char2):
        return match_award if char1 == char2 else sub_penalty

    if banded_width < 0:
        # Needleman-Wunsch implementation
        num_rows, num_cols = len(seq1) + 1, len(seq2) + 1
        val_table = [[0] * num_cols for _ in range(num_rows)]
        back_table = [[None] * num_cols for _ in range(num_rows)]

        # Initialize first row and column
        for i in range(1, num_rows):
            val_table[i][0] = i * indel_penalty
            back_table[i][0] = 'UP'
        for j in range(1, num_cols):
            val_table[0][j] = j * indel_penalty
            back_table[0][j] = 'LEFT'

        # Fill tables with tie-breaking
        for i in range(1, num_rows):
            for j in range(1, num_cols):
                match_score = val_table[i - 1][j - 1] + compare_chars(seq1[i - 1], seq2[j - 1])
                del_score = val_table[i - 1][j] + indel_penalty
                ins_score = val_table[i][j - 1] + indel_penalty

                # tie-breaking: DIAG > LEFT > UP
                if match_score <= ins_score and match_score <= del_score:
                    val_table[i][j] = match_score
                    back_table[i][j] = 'DIAG'
                elif ins_score < match_score and ins_score <= del_score:
                    val_table[i][j] = ins_score
                    back_table[i][j] = 'LEFT'
                else:
                    val_table[i][j] = del_score
                    back_table[i][j] = 'UP'

        # traceback
        alignment1, alignment2 = "", ""
        i, j = num_rows - 1, num_cols - 1
        while i > 0 or j > 0:
            if back_table[i][j] == 'DIAG':
                alignment1 = seq1[i - 1] + alignment1
                alignment2 = seq2[j - 1] + alignment2
                i -= 1
                j -= 1
            elif back_table[i][j] == 'LEFT':
                alignment1 = gap + alignment1
                alignment2 = seq2[j - 1] + alignment2
                j -= 1
            elif back_table[i][j] == 'UP':
                alignment1 = seq1[i - 1] + alignment1
                alignment2 = gap + alignment2
                i -= 1
        score = val_table[num_rows - 1][num_cols - 1]
        return score, alignment1, alignment2

    else:
        # Banded implementation
        n, m = len(seq1), len(seq2)
        if abs(n - m) > banded_width:
            return math.inf, None, None  

        band_size = 2 * banded_width + 1
        val_table = [[math.inf] * band_size for _ in range(n + 1)]
        back_table = [[None] * band_size for _ in range(n + 1)]

        # Initialize the diagonal of band
        val_table[0][banded_width] = 0
        back_table[0][banded_width] = 'START'

        # Initialize the first column of band
        for i in range(1, n + 1):
            k = banded_width - min(banded_width, i)
            val_table[i][k] = i * indel_penalty
            back_table[i][k] = 'UP'

        # Fill in just the band with tie-breaking
        for i in range(1, n + 1):
            for j in range(max(1, i - banded_width), min(m, i + banded_width) + 1):
                col = j - i + banded_width
                if col < 0 or col >= band_size:
                    continue
                diag_cost = val_table[i - 1][col] + compare_chars(seq1[i - 1], seq2[j - 1])
                left_cost = val_table[i][col - 1] + indel_penalty if col - 1 >= 0 else math.inf
                up_cost = val_table[i - 1][col + 1] + indel_penalty if col + 1 < band_size else math.inf

                # tie-breaking: DIAG > LEFT > UP
                if diag_cost <= left_cost and diag_cost <= up_cost:
                    val_table[i][col] = diag_cost
                    back_table[i][col] = 'DIAG'
                elif left_cost < diag_cost and left_cost <= up_cost:
                    val_table[i][col] = left_cost
                    back_table[i][col] = 'LEFT'
                else:
                    val_table[i][col] = up_cost
                    back_table[i][col] = 'UP'

        # Traceback
        alignment1, alignment2 = "", ""
        i, j = n, m
        col = j - i + banded_width
        while i > 0 or j > 0:
            if col < 0 or col >= band_size:
                break 

            direction = back_table[i][col]
            if direction == 'DIAG':
                alignment1 = seq1[i - 1] + alignment1
                alignment2 = seq2[j - 1] + alignment2
                i -= 1
                j -= 1
            elif direction == 'LEFT':
                alignment1 = gap + alignment1
                alignment2 = seq2[j - 1] + alignment2
                j -= 1
            elif direction == 'UP':
                alignment1 = seq1[i - 1] + alignment1
                alignment2 = gap + alignment2
                i -= 1
            else:
                break  

            col = j - i + banded_width  

        score = val_table[n][banded_width]
        return score, alignment1, alignment2
