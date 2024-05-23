void EXPECT_ARRAY_EQ(utils::four_vec x, utils::four_vec y, std::string msg ="")
{
    ASSERT_EQ(x.size(), y.size()) << "Vectors x and y are of unequal length" << msg;

    for (int i = 0; i < x.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(x[i], y[i]) << "Vectors x and y differ at index " << i << msg;
    }
}

void EXPECT_ARRAY_EQ(utils::r2_tensor x, utils::r2_tensor y, std::string msg ="")
{
    ASSERT_EQ(x.size(), y.size()) << "Vectors x and y are of unequal length" << msg;

    for (int i = 0; i < 4; ++i)
    {
        for (size_t j = 0; j < 4; j++)
        {
            EXPECT_DOUBLE_EQ(x[i][j], y[i][j]) << "Vectors x and y differ at index [" << i << ',' << j << "]" << msg;
        }
    }
}